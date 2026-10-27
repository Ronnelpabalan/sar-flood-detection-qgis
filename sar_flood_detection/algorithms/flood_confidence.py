"""Processing tool: flood confidence map from a pre-event and a post-event image, optionally with terrain."""

import numpy as np
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputNumber,
    QgsProcessingParameterBand,
    QgsProcessingParameterEnum,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
)

from ..core.confidence import HIGH, LOW, MEDIUM, NOT_FLOOD, PERMANENT_WATER, confidence_classes, flood_confidence
from ..core.thresholds import NODATA
from ..dem import warp_to_grid
from ..raster_io import cell_size_m, read_band, same_grid, write_bands, write_classes
from .common import (DEFAULT_WATER_THRESHOLD, FILTER_OPTIONS, UNIT_OPTIONS, driver_for, file_palette,
                     prepare_db, score_on_load, style_on_load)

CLASSES = [
    (NOT_FLOOD, "Not flood", (255, 255, 255, 0)),                     # transparent when loaded
    (LOW, "Flood, low confidence", (254, 217, 118, 255)),             # yellow
    (MEDIUM, "Flood, medium confidence", (253, 141, 60, 255)),        # orange
    (HIGH, "Flood, high confidence", (189, 0, 38, 255)),              # dark red
    (PERMANENT_WATER, "Permanent water", (31, 119, 180, 255)),        # blue
]

SCORE_COLORS = [
    (0.0, (255, 255, 255, 0), "0 (no flood evidence)"),
    (0.2, (254, 217, 118, 255), "0.2"),
    (0.5, (253, 141, 60, 255), "0.5 (old threshold)"),
    (0.8, (189, 0, 38, 255), "0.8"),
    (1.0, (103, 0, 13, 255), "1"),
]


class FloodConfidenceAlgorithm(QgsProcessingAlgorithm):
    PRE = "PRE"
    PRE_BAND = "PRE_BAND"
    POST = "POST"
    POST_BAND = "POST_BAND"
    UNITS = "UNITS"
    FILTER = "FILTER"
    WINDOW = "WINDOW"
    WATER_THRESHOLD = "WATER_THRESHOLD"
    MIN_DROP = "MIN_DROP"
    PRE_MASK = "PRE_MASK"
    POST_MASK = "POST_MASK"
    DEM = "DEM"
    MAX_SLOPE = "MAX_SLOPE"
    SLOPE_WINDOW = "SLOPE_WINDOW"
    NEAR_WATER = "NEAR_WATER"
    MAX_HEIGHT = "MAX_HEIGHT"
    OUTPUT = "OUTPUT"
    SCORE = "SCORE"
    HIGH_KM2 = "HIGH_KM2"
    MEDIUM_KM2 = "MEDIUM_KM2"
    LOW_KM2 = "LOW_KM2"
    PERMANENT_KM2 = "PERMANENT_KM2"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.PRE, "Pre-event raster"))
        self.addParameter(QgsProcessingParameterBand(self.PRE_BAND, "Pre-event band", 1, self.PRE))
        self.addParameter(QgsProcessingParameterRasterLayer(self.POST, "Post-event raster"))
        self.addParameter(QgsProcessingParameterBand(self.POST_BAND, "Post-event band", 1, self.POST))
        self.addParameter(QgsProcessingParameterEnum(
            self.UNITS, "Input units (both rasters)", UNIT_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum(
            self.FILTER, "Speckle filter (applied to both)", FILTER_OPTIONS, defaultValue=1))
        self.addParameter(QgsProcessingParameterNumber(
            self.WINDOW, "Filter window size (pixels, odd number)",
            QgsProcessingParameterNumber.Integer, defaultValue=5, minValue=3, maxValue=15))
        self.addParameter(QgsProcessingParameterNumber(
            self.WATER_THRESHOLD, "Water threshold (dB): score 0.5 here",
            QgsProcessingParameterNumber.Double, defaultValue=DEFAULT_WATER_THRESHOLD))
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_DROP, "Decrease for flooding (dB): score 0.5 here",
            QgsProcessingParameterNumber.Double, defaultValue=3.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.PRE_MASK, "Pre-event layover/shadow mask (optional)", optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.POST_MASK, "Post-event layover/shadow mask (optional)", optional=True))
        self.addParameter(QgsProcessingParameterMultipleLayers(
            self.DEM, "DEM for terrain checks (optional; one or more tiles, heights in m)",
            QgsProcessing.TypeRaster, optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.MAX_SLOPE, "Slope (degrees): score 0.5 here",
            QgsProcessingParameterNumber.Double, defaultValue=5.0, minValue=0.1))
        self.addParameter(QgsProcessingParameterNumber(
            self.SLOPE_WINDOW, "Distance over which the slope is taken (m)",
            QgsProcessingParameterNumber.Double, defaultValue=90.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.NEAR_WATER, "Ignore the slope up to this height above permanent water (m)",
            QgsProcessingParameterNumber.Double, defaultValue=3.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.MAX_HEIGHT, "Height above nearest permanent water (m): score 0.5 here",
            QgsProcessingParameterNumber.Double, defaultValue=10.0, minValue=0.1))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT, "Flood confidence classes"))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.SCORE, "Flood confidence score (0-1)", optional=True, createByDefault=True))
        for key, label in ((self.HIGH_KM2, "High confidence flood (km²)"),
                           (self.MEDIUM_KM2, "Medium confidence flood (km²)"),
                           (self.LOW_KM2, "Low confidence flood (km²)"),
                           (self.PERMANENT_KM2, "Permanent water (km²)")):
            self.addOutput(QgsProcessingOutputNumber(key, label))

    def _mask(self, parameters, name, context, grid):
        """Layover/shadow cells (any class above 0) of an optional mask, or None."""
        layer = self.parameterAsRasterLayer(parameters, name, context)
        if layer is None:
            return None
        mask, mask_grid = read_band(layer.source(), 1)
        if not same_grid(mask_grid, grid):
            raise QgsProcessingException(f"The {self.parameterDefinition(name).description()} is not on "
                                         "the same pixel grid as the images.")
        with np.errstate(invalid="ignore"):
            return mask > 0

    def processAlgorithm(self, parameters, context, feedback):
        pre_layer = self.parameterAsRasterLayer(parameters, self.PRE, context)
        post_layer = self.parameterAsRasterLayer(parameters, self.POST, context)
        if pre_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.PRE))
        if post_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.POST))
        units = self.parameterAsEnum(parameters, self.UNITS, context)
        filter_index = self.parameterAsEnum(parameters, self.FILTER, context)
        window = self.parameterAsInt(parameters, self.WINDOW, context)
        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)
        score_path = self.parameterAsOutputLayer(parameters, self.SCORE, context)

        pre, pre_grid = read_band(pre_layer.source(), self.parameterAsInt(parameters, self.PRE_BAND, context))
        post, grid = read_band(post_layer.source(), self.parameterAsInt(parameters, self.POST_BAND, context))
        if not same_grid(pre_grid, grid):
            raise QgsProcessingException(
                "The pre- and post-event rasters must have the same size, extent, pixel size and "
                "coordinate system. Align them first, e.g. with Raster > Align Rasters.")
        invalid = None
        for name in (self.PRE_MASK, self.POST_MASK):
            mask = self._mask(parameters, name, context, grid)
            if mask is not None:
                invalid = mask if invalid is None else invalid | mask
        feedback.setProgress(5)

        height = None
        dem_paths = [layer.source() for layer in self.parameterAsLayerList(parameters, self.DEM, context)]
        if dem_paths:
            feedback.pushInfo(f"Resampling {len(dem_paths)} DEM raster(s) onto the image grid")
            height = warp_to_grid(dem_paths, grid, fill=np.nan)
            covered = np.isfinite(height).mean()
            if covered < 0.99:
                feedback.pushWarning(f"The DEM covers only {100 * covered:.1f} % of the grid; "
                                     "terrain is not used where it is missing.")
        feedback.setProgress(15)

        pre = prepare_db(pre, units, filter_index, window, feedback, "to pre-event image")
        feedback.setProgress(40)
        post = prepare_db(post, units, filter_index, window, feedback, "to post-event image")
        feedback.setProgress(65)
        if feedback.isCanceled():
            return {}

        result = flood_confidence(
            pre, post, self.parameterAsDouble(parameters, self.WATER_THRESHOLD, context),
            self.parameterAsDouble(parameters, self.MIN_DROP, context), height=height,
            spacing=cell_size_m(grid), max_slope=self.parameterAsDouble(parameters, self.MAX_SLOPE, context),
            slope_window=self.parameterAsDouble(parameters, self.SLOPE_WINDOW, context),
            near_water=self.parameterAsDouble(parameters, self.NEAR_WATER, context),
            max_height=self.parameterAsDouble(parameters, self.MAX_HEIGHT, context), invalid=invalid)
        del pre, post, height
        classes = confidence_classes(result["score"], result["permanent"])
        feedback.setProgress(85)

        write_classes(out_path, classes, grid, file_palette(CLASSES), driver_for(out_path))
        style_on_load(context, out_path, CLASSES)
        outputs = {self.OUTPUT: out_path}
        if score_path:
            write_bands(score_path, [result["score"]], grid, names=["flood_confidence"])
            score_on_load(context, score_path, SCORE_COLORS)
            outputs[self.SCORE] = score_path
        feedback.setProgress(100)

        row_m, col_m = cell_size_m(grid)
        km2 = row_m * col_m / 1e6
        n_valid = max(int(np.sum(classes != NODATA)), 1)
        for key, value, label in ((self.HIGH_KM2, HIGH, "high confidence flood"),
                                  (self.MEDIUM_KM2, MEDIUM, "medium confidence flood"),
                                  (self.LOW_KM2, LOW, "low confidence flood"),
                                  (self.PERMANENT_KM2, PERMANENT_WATER, "permanent water")):
            n = int(np.sum(classes == value))
            outputs[key] = n * km2
            feedback.pushInfo(f"{label}: {n * km2:.1f} km² ({100.0 * n / n_valid:.2f} % of valid pixels)")
        return outputs

    def name(self):
        return "floodconfidence"

    def displayName(self):
        return "Flood confidence map (pre/post, optional terrain)"

    def group(self):
        return "Flood mapping"

    def groupId(self):
        return "floodmapping"

    def shortHelpString(self):
        return (
            "Scores each pixel from 0 to 1 for how strongly the evidence points to flooding, using a "
            "pre-event and a post-event Sentinel-1 image of the same area, track and polarisation.\n\n"
            "Evidence, each scored 0-1 and exactly 0.5 at its setting: the post-event image is dark "
            "(water threshold), the backscatter dropped (decrease for flooding), it was not water before; "
            "with a DEM also: the ground is flat (slope, taken over 90 m by default) and low (height "
            "above the nearest permanent water). The slope is ignored for cells within 3 m (default) of "
            "the water level: they can flood whatever the slope, and in a surface model such as "
            "Copernicus GLO-30 river banks, dikes and tree lines look steep. The score is the smallest of these, so the weakest evidence decides. Pixels in "
            "layover or shadow on either date (from the masks) get no score.\n\n"
            "Classes: 1 = low (score 0.2-0.5), 2 = medium (above 0.5), 3 = high (0.8 and above), "
            "4 = permanent water (dark before and after), 0 = not flood, 255 = no data. Without a DEM, "
            "medium + high is exactly the flood of the change-detection tool with the same settings.\n\n"
            "This is a confidence score, not a calibrated probability: it has not been fitted to "
            "reference flood maps. Any DEM in metres works (e.g. Copernicus GLO-30 tiles); only height "
            "differences and slopes are used, so the vertical datum does not matter."
        )

    def createInstance(self):
        return FloodConfidenceAlgorithm()
