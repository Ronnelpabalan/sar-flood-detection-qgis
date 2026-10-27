"""Shared options and steps for the Processing tools."""

import os

from qgis.core import (
    QgsPalettedRasterRenderer,
    QgsProcessingException,
    QgsProcessingLayerPostProcessorInterface,
    QgsRasterFileWriter,
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


# QGIS runs each tool on a temporary copy of the algorithm, so stylers are kept
# here (one per output path) until QGIS has loaded and styled the layer.
_STYLERS = {}


def style_on_load(context, path, classes):
    """Attach a ClassStyler to an output that QGIS will load when the tool finishes."""
    if context.willLoadLayerOnCompletion(path):
        _STYLERS[path] = ClassStyler(classes)
        context.layerToLoadOnCompletionDetails(path).setPostProcessor(_STYLERS[path])
