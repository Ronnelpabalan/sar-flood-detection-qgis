"""Range-Doppler terrain correction of a Sentinel-1 GRD product onto a map grid.

For every map cell: longitude, latitude and ellipsoidal height (DEM + geoid)
-> Range-Doppler model -> (line, pixel) in the radar image -> bilinear sample
of the calibrated backscatter. Layover and shadow are found on a coarser grid
(the DEM is 30 m anyway) by comparing points that fall in the same image lines.
"""

import numpy as np

from .core.geometry import geodetic_to_ecef
from .core.terrain import LAYOVER, SHADOW, bilinear, layover_shadow, mask_codes
from .core.thresholds import NODATA
from .dem import GridLonLat, ellipsoid_heights, make_grid

EARTH_RADIUS = 6371000.0                      # m, to turn ground angles into distances


def _inside(line, pixel, info):
    return (line >= 0) & (line <= info["lines"] - 1) & (pixel >= 0) & (pixel <= info["samples"] - 1)


class GeometryLattice:
    """Range-Doppler results for any grid cell and height, without solving each cell.

    The exact solution is computed on a lattice of grid nodes (every `step` cells)
    at three heights; values in between are bilinear in row/column and quadratic
    in height. The results change smoothly over hundreds of metres, so this is
    much faster than solving 40 million cells and stays close to exact (see tests).
    """

    def __init__(self, geo, lonlat, h_min, h_max, step=16, keys=("line", "pixel")):
        width, height = lonlat.grid["size"]
        self.rows = np.unique(np.r_[np.arange(0, height, step), height - 1]).astype("float64")
        self.cols = np.unique(np.r_[np.arange(0, width, step), width - 1]).astype("float64")
        if h_max - h_min < 2.0:
            h_min, h_max = h_min - 1.0, h_max + 1.0
        self.levels = np.array([h_min, (h_min + h_max) / 2, h_max])
        lon, lat = lonlat.at(self.rows, self.cols)
        tables = [geo.observe(geodetic_to_ecef(lat.ravel(), lon.ravel(), np.full(lon.size, h)))
                  for h in self.levels]
        self.tables = {k: np.stack([t[k].reshape(lon.shape) for t in tables]) for k in keys}

    @staticmethod
    def _weights(nodes, x):
        i = np.clip(np.searchsorted(nodes, x, side="right") - 1, 0, len(nodes) - 2)
        return i, (x - nodes[i]) / (nodes[i + 1] - nodes[i])

    def __call__(self, rows, cols, heights, keys=("line", "pixel")):
        """Dict of arrays shaped like heights (len(rows), len(cols)) for the given cells."""
        ir, fr = self._weights(self.rows, np.asarray(rows, "float64"))
        ic, fc = self._weights(self.cols, np.asarray(cols, "float64"))
        fr, fc = fr[:, None], fc[None, :]
        h0, h1, h2 = self.levels
        h = np.asarray(heights, "float64")
        lagrange = ((h - h1) * (h - h2) / ((h0 - h1) * (h0 - h2)),
                    (h - h0) * (h - h2) / ((h1 - h0) * (h1 - h2)),
                    (h - h0) * (h - h1) / ((h2 - h0) * (h2 - h1)))
        out = {}
        for key in keys:
            total = np.zeros(h.shape)
            for level, weight in enumerate(lagrange):
                t = self.tables[key][level]
                top = (1 - fc) * t[ir][:, ic] + fc * t[ir][:, ic + 1]
                bottom = (1 - fc) * t[ir + 1][:, ic] + fc * t[ir + 1][:, ic + 1]
                total += weight * ((1 - fr) * top + fr * bottom)
            out[key] = total
        return out


def _observe_blocks(geo, lonlat, heights, step, block_rows, feedback=None, lattice=None):
    """Range-Doppler view of every step-th grid cell, computed in blocks of rows to save memory.

    With a lattice, values are interpolated from it; otherwise each cell is solved exactly.
    """
    height, width = heights.shape
    keys = ("line", "pixel", "slant_range", "look_angle", "ground_angle")
    parts = []
    for r0 in range(0, height, block_rows * step):
        h = heights[r0:r0 + block_rows * step:step, ::step]
        if lattice is not None:
            rows = np.arange(r0, min(r0 + block_rows * step, height), step)
            v = lattice(rows, np.arange(0, width, step), h, keys)
        else:
            lon, lat = lonlat.block(r0, block_rows * step, col_step=step, row_step=step)
            v = geo.observe(geodetic_to_ecef(lat.ravel(), lon.ravel(), h.ravel()))
        parts.append({k: np.ravel(v[k]).astype("float64") for k in keys})
        if feedback is not None and feedback.isCanceled():
            raise RuntimeError("Canceled")
    return {k: np.concatenate([p[k] for p in parts]) for k in keys}


def coarse_pass(geo, info, lonlat, heights, spacing=10.0, mask_spacing=30.0, block_rows=256, feedback=None,
                lattice=None):
    """Layover/shadow mask (on the full grid) and the image window (col0, row0, width, height) to read.

    Works on every mask_spacing / spacing-th cell; the DEM is 30 m, so finer detail adds nothing.
    """
    height, width = heights.shape
    step = max(1, int(round(mask_spacing / spacing)))
    coarse_shape = heights[::step, ::step].shape
    v = _observe_blocks(geo, lonlat, heights, step, block_rows, feedback, lattice)
    inside = _inside(v["line"], v["pixel"], info)
    if not inside.any():
        raise ValueError("The area does not overlap this image.")
    lines_per_bin = max(1, int(round(mask_spacing / info["azimuth_spacing"])))
    flags = np.zeros(inside.size, np.uint8)
    flags[inside] = layover_shadow(np.floor(v["line"][inside] / lines_per_bin).astype(np.int64),
                                   v["ground_angle"][inside], v["look_angle"][inside],
                                   v["slant_range"][inside],
                                   separation=mask_spacing / EARTH_RADIUS)    # compare points >= one DEM cell apart
    coarse = mask_codes(flags, inside).reshape(coarse_shape)
    mask = np.repeat(np.repeat(coarse, step, 0), step, 1)[:height, :width]
    margin = 2 * step + 2
    col0 = max(int(np.floor(v["pixel"][inside].min())) - margin, 0)
    row0 = max(int(np.floor(v["line"][inside].min())) - margin, 0)
    col1 = min(int(np.ceil(v["pixel"][inside].max())) + margin + 1, info["samples"])
    row1 = min(int(np.ceil(v["line"][inside].max())) + margin + 1, info["lines"])
    return mask, (col0, row0, col1 - col0, row1 - row0)


def terrain_correct(product, pols, bbox, dem_paths, geoid_path, spacing=10.0, orbit_file=None,
                    mask_spacing=30.0, block_rows=256, remove_noise=True, exact=False, feedback=None):
    """Terrain-corrected sigma0 (linear) for each polarisation, plus a layover/shadow mask.

    product: SafeProduct; bbox: (lon_min, lat_min, lon_max, lat_max).
    exact: solve the Range-Doppler equations for every cell instead of using
    a GeometryLattice (slow; for checking).
    Returns {"grid", "bands": {pol: float32 array}, "mask": uint8 array, "window"}.
    """
    def progress(percent, text=None):
        if feedback is None:
            return
        if text:
            feedback.pushInfo(text)
        feedback.setProgress(percent)
        if feedback.isCanceled():
            raise RuntimeError("Canceled")

    grid = make_grid(bbox, spacing)
    width, height = grid["size"]
    progress(2, f"Map grid: {width} x {height} cells of {spacing:g} m, EPSG:{grid['epsg']}")
    lonlat = GridLonLat(grid)
    heights = ellipsoid_heights(dem_paths, geoid_path, grid)
    progress(8, f"Ellipsoidal heights: {np.nanmin(heights):.1f} to {np.nanmax(heights):.1f} m")

    info = product.info(pols[0])
    geo = product.geometry(pols[0], orbit_file)
    lattice = None
    if not exact:
        lattice = GeometryLattice(geo, lonlat, float(heights.min()), float(heights.max()),
                                  keys=("line", "pixel", "slant_range", "look_angle", "ground_angle"))

    # 1) coarse pass: where the area lies in the image, and layover/shadow
    mask, window = coarse_pass(geo, info, lonlat, heights, spacing, mask_spacing, block_rows, feedback, lattice)
    col0, row0 = window[0], window[1]
    col1, row1 = col0 + window[2], row0 + window[3]
    seen = mask[mask != NODATA]
    progress(20, f"Image window: columns {col0}-{col1}, lines {row0}-{row1}; "
                 f"layover {np.mean((seen & LAYOVER) > 0) * 100:.2f} %, "
                 f"shadow {np.mean((seen & SHADOW) > 0) * 100:.2f} % of the area")

    # 2) calibrated backscatter for that window
    images = {}
    for i, pol in enumerate(pols):
        images[pol] = product.backscatter(pol, window, remove_noise=remove_noise)
        progress(20 + 10 * (i + 1) // len(pols), f"Calibrated {pol}")

    # 3) fine pass: every map cell -> image position -> bilinear sample
    bands = {pol: np.full((height, width), np.nan, np.float32) for pol in pols}
    for r0 in range(0, height, block_rows):
        h = heights[r0:r0 + block_rows]
        if lattice is not None:
            v = lattice(np.arange(r0, r0 + h.shape[0]), np.arange(width), h)
            line, pixel = v["line"].ravel(), v["pixel"].ravel()
        else:
            lon, lat = lonlat.block(r0, block_rows)
            line, pixel = geo.locate(geodetic_to_ecef(lat.ravel(), lon.ravel(), h.ravel()))
        ok = _inside(line, pixel, info)
        rows, cols = np.where(ok, line - row0, np.nan), np.where(ok, pixel - col0, np.nan)
        for pol in pols:
            bands[pol][r0:r0 + h.shape[0]] = bilinear(images[pol], rows, cols).reshape(h.shape)
        progress(30 + int(70 * min(r0 + block_rows, height) / height))
    mask[~np.isfinite(bands[pols[0]])] = NODATA
    return {"grid": grid, "bands": bands, "mask": mask, "window": window}
