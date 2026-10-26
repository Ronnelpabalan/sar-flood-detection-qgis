"""Processing tool: water map from a single Sentinel-1 image."""

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingOutputNumber,
    QgsProcessingParameterBand,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
)

from ..core.thresholds import NODATA, classify_below, otsu_threshold
from ..raster_io import read_band, write_classes
from .common import DEFAULT_WATER_THRESHOLD, FILTER_OPTIONS, UNIT_OPTIONS, driver_for, prepare_db

WATER_COLORS = {0: (0, 0, 0, 0), 1: (31, 119, 180, 255)}      # not water transparent, water blue


class WaterMapAlgorithm(QgsProcessingAlgorithm):
    INPUT = "INPUT"
    BAND = "BAND"
    UNITS = "UNITS"
    FILTER = "FILTER"
    WINDOW = "WINDOW"
    METHOD = "METHOD"
    THRESHOLD = "THRESHOLD"
    OUTPUT = "OUTPUT"
    THRESHOLD_USED = "THRESHOLD_USED"
    WATER_PERCENT = "WATER_PERCENT"

    METHOD_OPTIONS = ["Fixed threshold", "Otsu (whole image)"]

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT, "Sentinel-1 backscatter raster"))
        self.addParameter(QgsProcessingParameterBand(
            self.BAND, "Band", 1, self.INPUT))
        self.addParameter(QgsProcessingParameterEnum(
            self.UNITS, "Input units", UNIT_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum(
            self.FILTER, "Speckle filter", FILTER_OPTIONS, defaultValue=1))
        self.addParameter(QgsProcessingParameterNumber(
            self.WINDOW, "Filter window size (pixels, odd number)",
            QgsProcessingParameterNumber.Integer, defaultValue=5, minValue=3, maxValue=15))
        self.addParameter(QgsProcessingParameterEnum(
            self.METHOD, "Threshold method", self.METHOD_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.THRESHOLD, "Fixed threshold (dB): water if backscatter is below",
            QgsProcessingParameterNumber.Double, defaultValue=DEFAULT_WATER_THRESHOLD))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.OUTPUT, "Water mask"))
        self.addOutput(QgsProcessingOutputNumber(self.THRESHOLD_USED, "Threshold used (dB)"))
        self.addOutput(QgsProcessingOutputNumber(self.WATER_PERCENT, "Water (% of valid pixels)"))

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT))
        band = self.parameterAsInt(parameters, self.BAND, context)
        units = self.parameterAsEnum(parameters, self.UNITS, context)
        filter_index = self.parameterAsEnum(parameters, self.FILTER, context)
        window = self.parameterAsInt(parameters, self.WINDOW, context)
        method = self.parameterAsEnum(parameters, self.METHOD, context)
        threshold = self.parameterAsDouble(parameters, self.THRESHOLD, context)
        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        arr, grid = read_band(layer.source(), band)
        feedback.setProgress(10)
        db = prepare_db(arr, units, filter_index, window, feedback)
        feedback.setProgress(60)
        if feedback.isCanceled():
            return {}

        if method == 1:
            otsu = otsu_threshold(db)
            threshold = otsu["threshold"]
            feedback.pushInfo(f"Otsu threshold: {threshold:.2f} dB "
                              f"(histogram valley ratio {otsu['valley_ratio']:.2f})")
            if not otsu["bimodal"]:
                feedback.pushWarning(
                    "The histogram does not show two clear peaks (water and land), so the Otsu "
                    "threshold is probably unreliable. Use a fixed threshold, or run Otsu on an "
                    "area that contains both water and land.")

        mask = classify_below(db, threshold)
        write_classes(out_path, mask, grid, WATER_COLORS, driver_for(out_path))
        feedback.setProgress(100)

        valid = mask != NODATA
        water_percent = 100.0 * float(mask[valid].mean()) if valid.any() else 0.0
        feedback.pushInfo(f"Threshold used: {threshold:.2f} dB. Water: {water_percent:.1f}% of valid pixels.")
        return {self.OUTPUT: out_path, self.THRESHOLD_USED: threshold, self.WATER_PERCENT: water_percent}

    def name(self):
        return "watermapsingleimage"

    def displayName(self):
        return "Water map from one image"

    def group(self):
        return "Flood mapping"

    def groupId(self):
        return "floodmapping"

    def shortHelpString(self):
        return (
            "Maps open water in a Sentinel-1 backscatter raster: optional speckle filter, "
            "then 'water if backscatter is below the threshold'. Output: 1 = water, "
            "0 = not water, 255 = no data.\n\n"
            "The default (median 5x5 filter, -23.3 dB) is the threshold that worked best for VH "
            "on the Sen1Floods11 training data. It is a starting point: check it against local "
            "data and adjust, especially for VV, which usually needs a higher value.\n\n"
            "Otsu chooses the threshold from the image itself. It needs both water and land in "
            "the image; a warning is shown when the histogram does not have two clear peaks.\n\n"
            "The whole band is processed in memory, so clip large scenes to the area of interest first."
        )

    def createInstance(self):
        return WaterMapAlgorithm()
