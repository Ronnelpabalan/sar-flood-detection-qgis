"""DEM tiles, geoid and the map grid for terrain correction.

DEM: Copernicus GLO-30 (30 m), public on AWS, one 1 x 1 degree tile per file:
https://copernicus-dem-30m.s3.amazonaws.com/<name>/<name>.tif
Heights are above the EGM2008 geoid; the radar geometry needs heights above the
WGS84 ellipsoid, so the geoid height (from the grid shipped with QGIS/PROJ) is added.
"""

import math
import os
import urllib.error
import urllib.request

import numpy as np
from osgeo import gdal, osr
from scipy.interpolate import RegularGridInterpolator

DEM_URL = "https://copernicus-dem-30m.s3.amazonaws.com"
GEOID_FILES = ("us_nga_egm08_25.tif", "us_nga_egm96_15.tif")     # EGM2008 first, EGM96 as fallback


def default_cache_dir():
    return os.path.join(os.path.expanduser("~"), ".sar_flood_detection", "dem")


def _fetch(url, timeout=300):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def tile_name(lat, lon):
    """GLO-30 tile whose south-west corner is (lat, lon), integers in degrees."""
    ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM"


def tile_names(bbox):
    """Tiles covering a lon/lat box (lon_min, lat_min, lon_max, lat_max)."""
    lon_min, lat_min, lon_max, lat_max = bbox
    lats = range(math.floor(lat_min), max(math.ceil(lat_max), math.floor(lat_min) + 1))
    lons = range(math.floor(lon_min), max(math.ceil(lon_max), math.floor(lon_min) + 1))
    return [tile_name(la, lo) for la in lats for lo in lons]


def find_dem_tiles(bbox, cache_dir=None, download=True, fetch=_fetch, feedback=None):
    """Local paths of the DEM tiles for a box, downloading missing ones into the cache.

    Tiles that do not exist on the server (open sea) are skipped; a small
    '.missing' file remembers them so they are not requested again.
    """
    cache_dir = cache_dir or default_cache_dir()
    paths, not_found = [], []
    for name in tile_names(bbox):
        path = os.path.join(cache_dir, name + ".tif")
        if os.path.exists(path):
            paths.append(path)
            continue
        if os.path.exists(path + ".missing"):
            continue
        if not download:
            not_found.append(name)
            continue
        if feedback is not None:
            feedback.pushInfo(f"Downloading DEM tile {name}")
        try:
            data = fetch(f"{DEM_URL}/{name}/{name}.tif")
        except urllib.error.HTTPError as e:
            if e.code not in (403, 404):
                raise
            os.makedirs(cache_dir, exist_ok=True)
            open(path + ".missing", "w").close()
            if feedback is not None:
                feedback.pushInfo(f"No DEM tile {name} (sea); heights there are set to 0 m")
            continue
        os.makedirs(cache_dir, exist_ok=True)
        with open(path + ".part", "wb") as f:
            f.write(data)
        os.replace(path + ".part", path)
        paths.append(path)
    if not_found:                     # never fall back to flat ground without saying so
        raise IOError(f"DEM tiles not found in {cache_dir} (downloading is off): " + ", ".join(not_found))
    return paths


def find_geoid():
    """Path of the EGM2008 (or EGM96) geoid grid in the PROJ data folders."""
    folders = list(osr.GetPROJSearchPaths() or [])
    for name in GEOID_FILES:
        for folder in folders:
            path = os.path.join(folder, name)
            if os.path.exists(path):
                return path
    raise IOError("No geoid grid found (us_nga_egm08_25.tif or us_nga_egm96_15.tif) in " + ", ".join(folders))


def utm_epsg(lon, lat):
    zone = min(int((lon + 180) // 6) + 1, 60)
    return (32600 if lat >= 0 else 32700) + zone


def _lonlat_to(epsg):
    src, dst = osr.SpatialReference(), osr.SpatialReference()
    src.ImportFromEPSG(4326)
    dst.ImportFromEPSG(epsg)
    for s in (src, dst):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)        # always x = lon, y = lat
    return osr.CoordinateTransformation(src, dst), osr.CoordinateTransformation(dst, src), dst


def make_grid(bbox, spacing=10.0, epsg=None):
    """Map grid covering a lon/lat box, in UTM by default.

    Edges are snapped to multiples of the spacing, so two images over the same
    box get exactly the same grid (needed for change detection).
    """
    lon_min, lat_min, lon_max, lat_max = bbox
    epsg = epsg or utm_epsg((lon_min + lon_max) / 2, (lat_min + lat_max) / 2)
    forward, _, srs = _lonlat_to(epsg)
    edge = np.linspace(0, 1, 21)
    lons = np.r_[lon_min + edge * (lon_max - lon_min), np.full(21, lon_max),
                 lon_max - edge * (lon_max - lon_min), np.full(21, lon_min)]
    lats = np.r_[np.full(21, lat_min), lat_min + edge * (lat_max - lat_min),
                 np.full(21, lat_max), lat_max - edge * (lat_max - lat_min)]
    xy = np.array(forward.TransformPoints(np.column_stack([lons, lats]).tolist()))
    x0 = math.floor(xy[:, 0].min() / spacing) * spacing
    x1 = math.ceil(xy[:, 0].max() / spacing) * spacing
    y0 = math.floor(xy[:, 1].min() / spacing) * spacing
    y1 = math.ceil(xy[:, 1].max() / spacing) * spacing
    return {
        "size": (int(round((x1 - x0) / spacing)), int(round((y1 - y0) / spacing))),
        "geotransform": (x0, spacing, 0.0, y1, 0.0, -spacing),
        "projection": srs.ExportToWkt(),
        "epsg": epsg,
    }


def warp_to_grid(paths, grid, fill=0.0, resampling="bilinear"):
    """Resample rasters (mosaicked) onto the grid as float32; cells with no data get fill."""
    width, height = grid["size"]
    x0, dx, _, y1, _, dy = grid["geotransform"]
    if not paths:
        return np.full((height, width), fill, np.float32)
    vrt = gdal.BuildVRT("/vsimem/sar_flood_mosaic.vrt", list(paths))
    out = gdal.Warp("", vrt, format="MEM", outputBounds=(x0, y1 + dy * height, x0 + dx * width, y1),
                    width=width, height=height, dstSRS=grid["projection"], resampleAlg=resampling,
                    outputType=gdal.GDT_Float32, dstNodata=np.nan, multithread=True)
    vrt = None
    gdal.Unlink("/vsimem/sar_flood_mosaic.vrt")
    data = out.GetRasterBand(1).ReadAsArray()
    out = None
    data[~np.isfinite(data)] = fill
    return data


def ellipsoid_heights(dem_paths, geoid_path, grid):
    """Heights above the WGS84 ellipsoid on the grid: DEM height (0 at sea) + geoid height."""
    heights = warp_to_grid(dem_paths, grid, fill=0.0)
    heights += warp_to_grid([geoid_path], grid, fill=0.0)
    return heights


class GridLonLat:
    """Longitude and latitude of grid cell centres.

    Exact transforms on a coarse lattice (every `step` cells), interpolated
    bilinearly in between; the map projection is smooth, so the error is about
    a millimetre at most (step 32 on a 10 m grid).
    """

    def __init__(self, grid, step=32):
        width, height = grid["size"]
        x0, dx, _, y1, _, dy = grid["geotransform"]
        self.grid = grid
        self.rows = np.unique(np.r_[np.arange(0, height, step), height - 1]).astype("float64")
        self.cols = np.unique(np.r_[np.arange(0, width, step), width - 1]).astype("float64")
        cc, rr = np.meshgrid(self.cols, self.rows)
        x, y = x0 + (cc + 0.5) * dx, y1 + (rr + 0.5) * dy
        _, inverse, _ = _lonlat_to(grid["epsg"])
        ll = np.array(inverse.TransformPoints(np.column_stack([x.ravel(), y.ravel()]).tolist()))
        self._lon = RegularGridInterpolator((self.rows, self.cols), ll[:, 0].reshape(x.shape))
        self._lat = RegularGridInterpolator((self.rows, self.cols), ll[:, 1].reshape(x.shape))
        self._lattice = ll

    def bbox(self, margin=0.01):
        """(lon_min, lat_min, lon_max, lat_max) of the whole grid plus a margin (degrees).

        A UTM rectangle reaches past the lon/lat box it was made from, so DEM
        tiles must be chosen from this, not from the box.
        """
        lon, lat = self._lattice[:, 0], self._lattice[:, 1]
        return lon.min() - margin, lat.min() - margin, lon.max() + margin, lat.max() + margin

    def block(self, row0, nrows, col_step=1, row_step=1):
        """(lon, lat) arrays for rows row0 .. row0 + nrows - 1 (every row_step), every col_step column."""
        width, height = self.grid["size"]
        rows = np.arange(row0, min(row0 + nrows, height), row_step, dtype="float64")
        cols = np.arange(0, width, col_step, dtype="float64")
        return self.at(rows, cols)

    def at(self, rows, cols):
        """(lon, lat) arrays of shape (len(rows), len(cols)) for the given row and column numbers."""
        cc, rr = np.meshgrid(np.asarray(cols, "float64"), np.asarray(rows, "float64"))
        points = np.column_stack([rr.ravel(), cc.ravel()])
        return self._lon(points).reshape(cc.shape), self._lat(points).reshape(cc.shape)
