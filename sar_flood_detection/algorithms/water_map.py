"""Processing tool: water map from a single Sentinel-1 image."""

from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBand,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
)


class WaterMapAlgorithm(QgsProcessingAlgorithm):
    INPUT = "INPUT"
    BAND = "BAND"
    UNITS = "UNITS"
    FILTER = "FILTER"
    WINDOW = "WINDOW"
    METHOD = "METHOD"
    THRESHOLD = "THRESHOLD"
    OUTPUT = "OUTPUT"

    UNIT_OPTIONS = ["Decibels (dB)", "Linear power"]
    FILTER_OPTIONS = ["None", "Median", "Lee"]
    METHOD_OPTIONS = ["Fixed threshold", "Otsu (whole image)"]

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT, "Sentinel-1 backscatter raster"))
        self.addParameter(QgsProcessingParameterBand(
            self.BAND, "Band", 1, self.INPUT))
        self.addParameter(QgsProcessingParameterEnum(
            self.UNITS, "Input units", self.UNIT_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum(
            self.FILTER, "Speckle filter", self.FILTER_OPTIONS, defaultValue=1))
        self.addParameter(QgsProcessingParameterNumber(
            self.WINDOW, "Filter window size (pixels, odd number)",
            QgsProcessingParameterNumber.Integer, defaultValue=5, minValue=3, maxValue=15))
        self.addParameter(QgsProcessingParameterEnum(
            self.METHOD, "Threshold method", self.METHOD_OPTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.THRESHOLD, "Fixed threshold (dB): water if backscatter is below",
            QgsProcessingParameterNumber.Double, defaultValue=-23.3))
        self.addParameter(QgsProcessingParameterRasterDestination(
            self.OUTPUT, "Water mask"))

    def processAlgorithm(self, parameters, context, feedback):
        raise QgsProcessingException("Not implemented yet.")

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
            "then 'water if backscatter is below the threshold'.\n\n"
            "The default (median 5x5 filter, -23.3 dB) is the threshold that worked best for VH "
            "on the Sen1Floods11 training data. It is a starting point: check it against local "
            "data and adjust, especially for VV, which usually needs a higher value."
        )

    def createInstance(self):
        return WaterMapAlgorithm()
