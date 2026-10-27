"""Shared options and steps for the Processing tools."""

import os

from qgis.core import (
    QgsColorRampShader,
    QgsContrastEnhancement,
    QgsPalettedRasterRenderer,
    QgsProcessingException,
    QgsProcessingLayerPostProcessorInterface,
    QgsRasterFileWriter,
    QgsRasterMinMaxOrigin,
    QgsRasterShader,
    QgsSingleBandGrayRenderer,
    QgsSingleBandPseudoColorRenderer,
)
from qgis.PyQt.QtGui import QColor

from ..core.filters import lee_filter_db, median_filter_db
from ..core.thresholds import to_db

UNIT_OPTIONS = ["Decibels (dB)", "Linear power"]
FILTER_OPTIONS = ["None", "Median", "Lee"]
DEFAULT_WATER_THRESHOLD = -23.3     # VH, median 5x5, best on Sen1Floods11 train (sar-flood-unet)


def prepare_db(arr, units, filter_index, window, feedback=None, label=""):
    """Convert to dB if needed and apply the chosen speckle filter."""
    db = to_db(arr) if units == 1 else arr
    if filter_index == 0:
        return db
    if window % 2 == 0:
        raise QgsProcessingException(f"Filter window size must be an odd number, got {window}.")
    if feedback is not None:
        feedback.pushInfo(f"Applying {FILTER_OPTIONS[filter_index]} filter {window}x{window} {label}".strip())
    if filter_index == 1:
        return median_filter_db(db, size=window)
    return lee_filter_db(db, size=window)


def driver_for(path):
    """GDAL driver matching the output file extension (GeoTIFF if unknown)."""
    return QgsRasterFileWriter.driverForExtension(os.path.splitext(path)[1]) or "GTiff"


def file_palette(classes):
    """Colour table for the output file: {value: (r, g, b, 255)}.

    GeoTIFF palettes cannot store transparency, so transparent classes are
    written as opaque colours here and made transparent by ClassStyler.
    """
    return {value: (*rgba[:3], 255) for value, _, rgba in classes}


class ClassStyler(QgsProcessingLayerPostProcessorInterface):
    """Styles a class raster with named, possibly transparent classes when QGIS loads it.

    classes: list of (value, label, (red, green, blue, alpha))
    """

    def __init__(self, classes):
        super().__init__()
        self.classes = classes

    def postProcessLayer(self, layer, context, feedback):
        items = [QgsPalettedRasterRenderer.Class(value, QColor(*rgba), label)
                 for value, label, rgba in self.classes]
        layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, items))
        layer.triggerRepaint()


class GrayStyler(QgsProcessingLayerPostProcessorInterface):
    """Shows one band of a raster in grey, stretched between its 2 % and 98 % values."""

    def __init__(self, band, lower=0.02, upper=0.98):
        super().__init__()
        self.band, self.lower, self.upper = band, lower, upper

    def postProcessLayer(self, layer, context, feedback):
        renderer = QgsSingleBandGrayRenderer(layer.dataProvider(), self.band)
        origin = renderer.minMaxOrigin()
        origin.setLimits(QgsRasterMinMaxOrigin.CumulativeCut)
        origin.setCumulativeCutLower(self.lower)
        origin.setCumulativeCutUpper(self.upper)
        renderer.setMinMaxOrigin(origin)
        layer.setRenderer(renderer)
        layer.setContrastEnhancement(QgsContrastEnhancement.StretchToMinimumMaximum,
                                     QgsRasterMinMaxOrigin.CumulativeCut)
        layer.triggerRepaint()


class ScoreStyler(QgsProcessingLayerPostProcessorInterface):
    """Shows a 0-1 score with a colour ramp; colours: [(value, (red, green, blue, alpha), label)]."""

    def __init__(self, colors):
        super().__init__()
        self.colors = colors

    def postProcessLayer(self, layer, context, feedback):
        ramp = QgsColorRampShader(self.colors[0][0], self.colors[-1][0])
        ramp.setColorRampType(QgsColorRampShader.Interpolated)
        ramp.setColorRampItemList([QgsColorRampShader.ColorRampItem(value, QColor(*rgba), label)
                                   for value, rgba, label in self.colors])
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(ramp)
        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)
        renderer.setClassificationMin(self.colors[0][0])
        renderer.setClassificationMax(self.colors[-1][0])
        layer.setRenderer(renderer)
        layer.triggerRepaint()


# QGIS runs each tool on a temporary copy of the algorithm, so stylers are kept
# here (one per output path) until QGIS has loaded and styled the layer.
_STYLERS = {}


def _attach(context, path, styler):
    if context.willLoadLayerOnCompletion(path):
        _STYLERS[path] = styler
        context.layerToLoadOnCompletionDetails(path).setPostProcessor(styler)


def style_on_load(context, path, classes):
    """Attach a ClassStyler to an output that QGIS will load when the tool finishes."""
    _attach(context, path, ClassStyler(classes))


def gray_on_load(context, path, band):
    """Show one band of an output in grey (2-98 % stretch) when QGIS loads it."""
    _attach(context, path, GrayStyler(band))


def score_on_load(context, path, colors):
    """Show a 0-1 score with a colour ramp when QGIS loads it."""
    _attach(context, path, ScoreStyler(colors))
