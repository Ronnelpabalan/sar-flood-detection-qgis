"""Reading and writing rasters with GDAL (included with QGIS)."""

import numpy as np
from osgeo import gdal, osr

from .core.thresholds import NODATA


def read_band(path, band=1):
    """Read one band as float32 with nodata values set to NaN.

    Returns (array, grid), where grid holds the size, geotransform and
    projection needed to write an output on the same pixel grid.
    """
    ds = gdal.Open(path)
    if ds is None:
        raise IOError(f"Cannot open raster: {path}")
    if not 1 <= band <= ds.RasterCount:
        raise ValueError(f"Band {band} requested, but the raster has {ds.RasterCount} band(s).")
    b = ds.GetRasterBand(band)
    arr = b.ReadAsArray().astype("float32")
    nodata = b.GetNoDataValue()
    if nodata is not None and not np.isnan(nodata):
        arr[arr == nodata] = np.nan
    grid = {
        "size": (ds.RasterXSize, ds.RasterYSize),
        "geotransform": ds.GetGeoTransform(),
        "projection": ds.GetProjection(),
    }
    return arr, grid


def same_grid(a, b):
    """True if two grids have the same size, pixel positions and coordinate system."""
    if a["size"] != b["size"]:
        return False
    if not np.allclose(a["geotransform"], b["geotransform"], rtol=1e-9, atol=1e-9):
        return False
    sa, sb = osr.SpatialReference(), osr.SpatialReference()
    sa.ImportFromWkt(a["projection"])
    sb.ImportFromWkt(b["projection"])
    return bool(sa.IsSame(sb))


def write_bands(path, bands, grid, names=None):
    """Write float32 bands (list of 2-D arrays) with NaN as nodata; names become band descriptions."""
    width, height = grid["size"]
    ds = gdal.GetDriverByName("GTiff").Create(path, width, height, len(bands), gdal.GDT_Float32,
                                              options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER"])
    if ds is None:
        raise IOError(f"Cannot create raster: {path}")
    ds.SetGeoTransform(grid["geotransform"])
    ds.SetProjection(grid["projection"])
    for i, data in enumerate(bands):
        b = ds.GetRasterBand(i + 1)
        b.SetNoDataValue(float("nan"))
        if names:
            b.SetDescription(names[i])
        b.WriteArray(data)
    ds.FlushCache()
    ds = None


def write_classes(path, classes, grid, colors, driver="GTiff"):
    """Write a uint8 class raster with NODATA set and a colour table.

    colors: {class value: (red, green, blue, alpha)}
    """
    width, height = grid["size"]
    options = ["COMPRESS=DEFLATE", "TILED=YES"] if driver == "GTiff" else []
    ds = gdal.GetDriverByName(driver).Create(path, width, height, 1, gdal.GDT_Byte, options=options)
    if ds is None:
        raise IOError(f"Cannot create raster: {path}")
    ds.SetGeoTransform(grid["geotransform"])
    ds.SetProjection(grid["projection"])
    b = ds.GetRasterBand(1)
    table = gdal.ColorTable()                       # GeoTIFF needs the palette before the pixel data
    for value, rgba in colors.items():
        table.SetColorEntry(value, rgba)
    b.SetRasterColorTable(table)
    b.SetNoDataValue(NODATA)
    b.WriteArray(classes)
    ds.FlushCache()
    ds = None
