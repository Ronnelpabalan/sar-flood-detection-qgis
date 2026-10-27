"""Processing tool: terrain-corrected backscatter from a raw Sentinel-1 GRD product (.SAFE.zip)."""

import os
import urllib.error

import numpy as np
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterEnum,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFile,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
)

from ..core.orbit import parse_time
from ..core.terrain import LAYOVER, SHADOW
from ..core.thresholds import NODATA
from ..dem import GridLonLat, find_dem_tiles, find_geoid, make_grid
from ..orbit_files import find_orbit_file
from ..raster_io import write_bands, write_classes
from ..safe import SafeProduct
from ..terrain_correction import terrain_correct
from .common import file_palette, gray_on_load, style_on_load

POL_OPTIONS = ["VV and VH", "VV only", "VH only"]
POL_CHOICES = [["VV", "VH"], ["VV"], ["VH"]]
ORBIT_OPTIONS = ["Precise orbit file (downloaded from ESA if needed)", "Orbit stored in the product"]
OUTPUT_UNITS = ["Decibels (dB)", "Linear power"]
MASK_CLASSES = [
    (0, "Seen normally", (255, 255, 255, 0)),                    # transparent when loaded
    (LAYOVER, "Layover", (255, 140, 0, 255)),                     # orange
    (SHADOW, "Shadow", (90, 90, 90, 255)),                        # grey
    (LAYOVER | SHADOW, "Layover and shadow", (200, 0, 200, 255)), # purple
]


class PrepareGrdAlgorithm(QgsProcessingAlgorithm):
    INPUT = "INPUT"
    EXTENT = "EXTENT"
    POLARISATIONS = "POLARISATIONS"
    SPACING = "SPACING"
    ORBIT = "ORBIT"
    DATA_FOLDER = "DATA_FOLDER"
    DOWNLOAD = "DOWNLOAD"
    REMOVE_NOISE = "REMOVE_NOISE"
    UNITS = "UNITS"
    OUTPUT = "OUTPUT"
    MASK = "MASK"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(
            self.INPUT, "Sentinel-1 GRD product (.SAFE.zip, or manifest.safe of an unzipped product)",
            fileFilter="Sentinel-1 products (*.zip *.ZIP *.safe *.SAFE)"))
        self.addParameter(QgsProcessingParameterExtent(self.EXTENT, "Area to process"))
        self.addParameter(QgsProcessingParameterEnum(
            self.POLARISATIONS, "Polarisations", POL_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.SPACING, "Pixel size of the output (m)", QgsProcessingParameterNumber.Double,
            defaultValue=10.0, minValue=10.0, maxValue=1000.0))
        self.addParameter(QgsProcessingParameterEnum(self.ORBIT, "Orbit", ORBIT_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterFile(
            self.DATA_FOLDER, "Folder for DEM tiles and orbit files (empty = ~/.sar_flood_detection)",
            behavior=QgsProcessingParameterFile.Folder, optional=True))
        self.addParameter(QgsProcessingParameterBoolean(
            self.DOWNLOAD, "Download missing DEM tiles and orbit files", defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(
            self.REMOVE_NOISE, "Remove thermal noise", defaultValue=True))
        self.addParameter(QgsProcessingParameterEnum(
            self.UNITS, "Output units", OUTPUT_UNITS, defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.OUTPUT, "Terrain-corrected backscatter (sigma0)"))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.MASK, "Layover/shadow mask", optional=True, createByDefault=True))

    def processAlgorithm(self, parameters, context, feedback):
        path = self.parameterAsFile(parameters, self.INPUT, context)
        try:
            product = SafeProduct(path)
        except (IOError, OSError) as e:
            raise QgsProcessingException(str(e))
        pols = POL_CHOICES[self.parameterAsEnum(parameters, self.POLARISATIONS, context)]
        missing = [p for p in pols if p not in product.polarisations]
        if missing:
            raise QgsProcessingException(
                f"The product has no {', '.join(missing)} (it has {', '.join(product.polarisations)}).")
        extent = self.parameterAsExtent(parameters, self.EXTENT, context, QgsCoordinateReferenceSystem("EPSG:4326"))
        bbox = (extent.xMinimum(), extent.yMinimum(), extent.xMaximum(), extent.yMaximum())
        if extent.isEmpty():
            raise QgsProcessingException("The area to process is empty.")
        spacing = self.parameterAsDouble(parameters, self.SPACING, context)
        use_orbit_file = self.parameterAsEnum(parameters, self.ORBIT, context) == 0
        folder = self.parameterAsFile(parameters, self.DATA_FOLDER, context) or None
        download = self.parameterAsBool(parameters, self.DOWNLOAD, context)
        remove_noise = self.parameterAsBool(parameters, self.REMOVE_NOISE, context)
        in_db = self.parameterAsEnum(parameters, self.UNITS, context) == 0
        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)
        mask_path = self.parameterAsOutputLayer(parameters, self.MASK, context)

        info = product.info(pols[0])
        feedback.pushInfo(f"{info['mission']} {info['pass'].lower()} pass, {info['start_time']}, "
                          f"polarisations {', '.join(pols)}")
        feedback.pushInfo("Area: lon {:.4f} to {:.4f}, lat {:.4f} to {:.4f}".format(bbox[0], bbox[2], bbox[1], bbox[3]))

        orbit_file = None
        if use_orbit_file:
            orbit_dir = os.path.join(folder, "orbits") if folder else None
            try:
                orbit_file = find_orbit_file(info["mission"], parse_time(info["start_time"]),
                                             parse_time(info["stop_time"]), orbit_dir, download, feedback=feedback)
            except (urllib.error.URLError, OSError) as e:
                feedback.pushWarning(f"Could not get an orbit file ({e}).")
            if orbit_file:
                feedback.pushInfo(f"Orbit file: {os.path.basename(orbit_file)}")
            else:
                feedback.pushWarning("No orbit file found; using the orbit stored in the product.")

        grid = make_grid(bbox, spacing)
        dem_dir = os.path.join(folder, "dem") if folder else None
        try:
            dem_paths = find_dem_tiles(GridLonLat(grid).bbox(), dem_dir, download, feedback=feedback)
            geoid = find_geoid()
        except (urllib.error.URLError, OSError) as e:
            raise QgsProcessingException(f"Cannot get the DEM or geoid: {e}")
        feedback.pushInfo(f"DEM tiles: {len(dem_paths)}; geoid: {os.path.basename(geoid)}")
        if feedback.isCanceled():
            return {}

        try:
            result = terrain_correct(product, pols, bbox, dem_paths, geoid, spacing, orbit_file,
                                     remove_noise=remove_noise, feedback=feedback)
        except ValueError as e:
            raise QgsProcessingException(str(e))
        except RuntimeError:                                  # canceled by the user
            return {}

        bands, names = [], []
        for pol in pols:
            band = result["bands"][pol]
            if in_db:
                with np.errstate(divide="ignore", invalid="ignore"):
                    np.log10(band, out=band)
                band *= 10.0
            bands.append(band)
            names.append(f"Sigma0_{pol}" + ("_dB" if in_db else ""))
        write_bands(out_path, bands, result["grid"], names)
        gray_on_load(context, out_path, pols.index("VH") + 1 if "VH" in pols else 1)   # VH shows water best
        outputs = {self.OUTPUT: out_path}
        if mask_path:
            write_classes(mask_path, result["mask"], result["grid"], file_palette(MASK_CLASSES))
            style_on_load(context, mask_path, MASK_CLASSES)
            outputs[self.MASK] = mask_path

        mask = result["mask"]
        seen = mask[mask != NODATA]
        if seen.size:
            feedback.pushInfo(f"Valid cells: {100.0 * seen.size / mask.size:.1f} %; layover "
                              f"{100.0 * np.mean((seen & LAYOVER) > 0):.2f} %, shadow "
                              f"{100.0 * np.mean((seen & SHADOW) > 0):.2f} %.")
        feedback.pushInfo("Bands: " + ", ".join(f"{i + 1} = {n}" for i, n in enumerate(names)))
        return outputs

    def name(self):
        return "preparegrd"

    def displayName(self):
        return "Prepare Sentinel-1 GRD (from SAFE zip)"

    def group(self):
        return "Preprocessing"

    def groupId(self):
        return "preprocessing"

    def shortHelpString(self):
        return (
            "Turns a raw Sentinel-1 GRD product (the .SAFE.zip as downloaded) into terrain-corrected "
            "backscatter on a map grid, without other software.\n\n"
            "Steps: calibration to sigma0 and thermal noise removal; precise orbit (from ESA, cached); "
            "Copernicus GLO-30 DEM tiles (from AWS, cached) plus the EGM2008 geoid shipped with QGIS; "
            "Range-Doppler terrain correction onto a UTM grid; bilinear resampling.\n\n"
            "Outputs: one band per polarisation (dB or linear), NaN = no data; and a layover/shadow "
            "mask: 1 = layover, 2 = shadow, 3 = both, 0 = seen normally, 255 = no data. Results in "
            "layover and shadow areas are not reliable.\n\n"
            "For change detection, process both dates with exactly the same area and pixel size: "
            "the grid edges are rounded to whole pixels, so the outputs then line up cell for cell.\n\n"
            "DEM tiles and orbit files are kept in the chosen folder (sub-folders dem and orbits), so "
            "they are downloaded only once. A 60 x 65 km area at 10 m takes about 2 minutes and "
            "1.3 GB of memory."
        )

    def createInstance(self):
        return PrepareGrdAlgorithm()
