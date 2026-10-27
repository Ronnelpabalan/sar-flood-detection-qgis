"""Terrain-correction helpers: resampling the radar image and finding layover and shadow.

Mask values: 0 = seen normally, 1 = layover, 2 = shadow, 3 = both, NODATA = outside the image.
"""

import numpy as np

from .thresholds import NODATA

LAYOVER, SHADOW = 1, 2


def bilinear(image, rows, cols):
    """Bilinear sample of image at fractional (rows, cols); NaN outside or next to NaN pixels."""
    rows, cols = np.asarray(rows, dtype="float64"), np.asarray(cols, dtype="float64")
    h, w = image.shape
    ok = np.isfinite(rows) & np.isfinite(cols) & (rows >= 0) & (cols >= 0) & (rows <= h - 1) & (cols <= w - 1)
    r = np.where(ok, rows, 0.0)
    c = np.where(ok, cols, 0.0)
    r0 = np.minimum(np.floor(r).astype(np.int64), max(h - 2, 0))
    c0 = np.minimum(np.floor(c).astype(np.int64), max(w - 2, 0))
    fr, fc = r - r0, c - c0
    r1, c1 = np.minimum(r0 + 1, h - 1), np.minimum(c0 + 1, w - 1)
    out = ((1 - fr) * (1 - fc) * image[r0, c0] + (1 - fr) * fc * image[r0, c1]
           + fr * (1 - fc) * image[r1, c0] + fr * fc * image[r1, c1])
    out = out.astype(image.dtype, copy=False)
    out[~ok] = np.nan
    return out


def _nearer_max(values, rank, position, separation):
    """Largest value among the points of the same group that are at least
    `separation` nearer to the track (-inf if there are none).

    Input is sorted by rank (group number, non-decreasing), then position.
    Adding rank * span lifts each group above the previous one, so one
    cumulative maximum restarts at every group.
    """
    v = values.astype("float64")
    low = v.min()
    span = v.max() - low + 1.0
    best = np.maximum.accumulate(v - low + rank * span)
    p_span = position.max() - position.min() + separation + 1.0
    key = position - position.min() + rank * p_span
    j = np.searchsorted(key, key - separation, side="right") - 1   # last point far enough behind
    j_safe = np.maximum(j, 0)
    ok = (j >= 0) & (rank[j_safe] == rank)
    return np.where(ok, best[j_safe] - rank * span + low, -np.inf)


def _farther_min(values, rank, position, separation):
    """Smallest value among the points of the same group at least `separation` farther away."""
    flip = slice(None, None, -1)
    return -_nearer_max(-values[flip], rank.max() - rank[flip], -position[flip], separation)[flip]


def layover_shadow(group, ground_angle, look_angle, slant_range, separation=0.0, tolerance=1e-9):
    """Layover/shadow flags for points seen in the same image lines.

    group: image line bin of each point. Within a bin, points are ordered by
    distance from the ground track (ground_angle). A point is in shadow if a
    nearer point has a larger look angle (the beam hits that point first). Two
    points are in layover if the farther one is closer to the satellite
    (smaller slant range), so their echoes arrive mixed together.

    Only points at least `separation` (same unit as ground_angle) apart are
    compared. Without it, two neighbours at almost the same distance would be
    flagged for any small height difference, e.g. trees in a surface model.
    Returns uint8 flags (LAYOVER | SHADOW) in the input order.
    """
    group = np.asarray(group)
    flags = np.zeros(len(group), np.uint8)
    if len(group) == 0:
        return flags
    order = np.lexsort((ground_angle, group))
    g = group[order]
    rank = np.cumsum(np.r_[0, np.diff(g) != 0]).astype("float64")
    look, rng, pos = look_angle[order], slant_range[order], np.asarray(ground_angle, "float64")[order]
    flags[look < _nearer_max(look, rank, pos, separation) - tolerance] |= SHADOW
    layover = ((rng < _nearer_max(rng, rank, pos, separation) - 1e-3)
               | (rng > _farther_min(rng, rank, pos, separation) + 1e-3))
    flags[layover] |= LAYOVER
    out = np.empty_like(flags)
    out[order] = flags
    return out


def mask_codes(flags, valid):
    """uint8 mask from flags, with NODATA where the point is outside the image."""
    out = flags.astype(np.uint8, copy=True)
    out[~valid] = NODATA
    return out
