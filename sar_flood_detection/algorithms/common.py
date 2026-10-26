"""Shared options and steps for the Processing tools."""

import os

from qgis.core import QgsProcessingException, QgsRasterFileWriter

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
