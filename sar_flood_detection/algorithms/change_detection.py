"""Processing tool: flood map from a pre-event and a post-event Sentinel-1 image."""

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

from ..core.change import DRY, FLOOD, PERMANENT_WATER, change_detection
from ..core.thresholds import NODATA
from ..raster_io import read_band, same_grid, write_classes
from .common import DEFAULT_WATER_THRESHOLD, FILTER_OPTIONS, UNIT_OPTIONS, driver_for, prepare_db

CLASS_COLORS = {
    DRY: (0, 0, 0, 0),                      # transparent
    FLOOD: (214, 39, 40, 255),              # red
    PERMANENT_WATER: (31, 119, 180, 255),   # blue
}


class ChangeDetectionAlgorithm(QgsProcessingAlgorithm):
    PRE = "PRE"
    PRE_BAND = "PRE_BAND"
    POST = "POST"
    POST_BAND = "POST_BAND"
    UNITS = "UNITS"
    FILTER = "FILTER"
    WINDOW = "WINDOW"
    WATER_THRESHOLD = "WATER_THRESHOLD"
    MIN_DROP = "MIN_DROP"
    OUTPUT = "OUTPUT"
    FLOOD_PERCENT = "FLOOD_PERCENT"
    PERMANENT_PERCENT = "PERMANENT_PERCENT"

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
            self.WATER_THRESHOLD, "Water threshold (dB): water if backscatter is below",
            QgsProcessingParameterNumber.Double, defaultValue=DEFAULT_WATER_THRESHOLD))
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_DROP, "Minimum decrease for flooding (dB)",
            QgsProcessingParameterNumber.Double, defaultValue=3.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT, "Flood map"))
        self.addOutput(QgsProcessingOutputNumber(self.FLOOD_PERCENT, "Flood (% of valid pixels)"))
        self.addOutput(QgsProcessingOutputNumber(self.PERMANENT_PERCENT, "Permanent water (% of valid pixels)"))

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
        water_threshold = self.parameterAsDouble(parameters, self.WATER_THRESHOLD, context)
        min_drop = self.parameterAsDouble(parameters, self.MIN_DROP, context)
        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        pre, pre_grid = read_band(pre_layer.source(), self.parameterAsInt(parameters, self.PRE_BAND, context))
        post, post_grid = read_band(post_layer.source(), self.parameterAsInt(parameters, self.POST_BAND, context))
        if not same_grid(pre_grid, post_grid):
            raise QgsProcessingException(
                "The pre- and post-event rasters must have the same size, extent, pixel size and "
                "coordinate system. Align them first, e.g. with Raster > Align Rasters.")
        feedback.setProgress(10)

        pre_db = prepare_db(pre, units, filter_index, window, feedback, "to pre-event image")
        feedback.setProgress(40)
        post_db = prepare_db(post, units, filter_index, window, feedback, "to post-event image")
        feedback.setProgress(70)
        if feedback.isCanceled():
            return {}

        classes = change_detection(pre_db, post_db, water_threshold, min_drop)
        write_classes(out_path, classes, post_grid, CLASS_COLORS, driver_for(out_path))
        feedback.setProgress(100)

        valid = classes != NODATA
        n_valid = max(int(valid.sum()), 1)
        flood_percent = 100.0 * int((classes == FLOOD).sum()) / n_valid
        permanent_percent = 100.0 * int((classes == PERMANENT_WATER).sum()) / n_valid
        feedback.pushInfo(f"Flood: {flood_percent:.1f}% and permanent water: {permanent_percent:.1f}% "
                          "of valid pixels.")
        return {self.OUTPUT: out_path, self.FLOOD_PERCENT: flood_percent,
                self.PERMANENT_PERCENT: permanent_percent}

    def name(self):
        return "floodchangedetection"

    def displayName(self):
        return "Flood map from change detection (pre/post)"

    def group(self):
        return "Flood mapping"

    def groupId(self):
        return "floodmapping"

    def shortHelpString(self):
        return (
            "Compares a pre-event and a post-event Sentinel-1 image of the same area and polarisation.\n\n"
            "Output classes: 1 = flood (dark after the event, not before, and backscatter dropped by "
            "more than the minimum decrease), 2 = permanent water (dark before and after), "
            "0 = dry or no clear change, 255 = no data.\n\n"
            "Both rasters must be on the same pixel grid (same size, extent and coordinate system). "
            "Use images from the same orbit track and polarisation, so that differences come from "
            "the flood rather than from the viewing geometry.\n\n"
            "The default water threshold (-23.3 dB) suits VH after a median 5x5 filter; VV usually "
            "needs a higher value. Flooded vegetation can become brighter rather than darker and is "
            "not detected by this rule."
        )

    def createInstance(self):
        return ChangeDetectionAlgorithm()
