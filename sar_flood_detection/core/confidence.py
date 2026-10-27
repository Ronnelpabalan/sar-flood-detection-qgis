"""Flood confidence score from a pre-event and a post-event image, optionally with terrain.

Each piece of evidence becomes a score from 0 to 1 that rises smoothly
(a logistic curve) and is exactly 0.5 at the usual threshold:

    1  the post-event image is dark       0.5 at post = water_threshold
    2  the backscatter dropped            0.5 at a drop of min_drop
    3  it was not water before the event  0.5 at pre = water_threshold
    5a the ground is flat (optional)      0.5 at max_slope degrees, with the slope taken over
                                          slope_window metres; ignored for cells less than
                                          near_water metres above the nearest permanent water
    5b the ground is low (optional)       0.5 at max_height metres above the nearest permanent water

(4, layover and shadow, is handled by the invalid mask: those pixels get no score.)

The flood score is the smallest of these scores: every condition must hold,
so the weakest evidence decides. Without terrain, a score above 0.5 means
the same pixels as change_detection() (apart from values exactly on a
threshold); terrain can only lower the score.

Why the slope is ignored near the water: a cell at the water level can
flood whatever the slope there. The Copernicus DEM is a surface model, so
river banks, dikes and tree lines look steep; on Thanh Hoa 2025 a plain
10 m slope test mostly removed dark cells on river banks (median 50 m from
the river, 1.8 m above it).

This is a confidence score, not a calibrated probability: it has not been
fitted to reference flood maps.

Classes: NOT_FLOOD = 0, LOW = 1, MEDIUM = 2, HIGH = 3, PERMANENT_WATER = 4, NODATA = 255.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt, uniform_filter

from .thresholds import NODATA

NOT_FLOOD = 0
LOW = 1
MEDIUM = 2
HIGH = 3
PERMANENT_WATER = 4


def rising(x, centre, width):
    """Logistic curve: 0.5 at x = centre, about 0.88 at centre + 2 * width. NaN stays NaN."""
    z = (np.asarray(x, dtype="float32") - np.float32(centre)) / np.float32(width)
    return (1.0 / (1.0 + np.exp(-np.clip(z, -50.0, 50.0)))).astype("float32")


def slope_degrees(height, spacing, window=None):
    """Terrain slope in degrees; spacing is the cell size in metres, or a (row, column) pair.

    window: optional distance in metres over which heights are averaged first, so the slope
    describes the land form rather than single trees, houses or dikes in a surface model.
    """
    row_m, col_m = (spacing, spacing) if np.isscalar(spacing) else spacing
    height = np.asarray(height, dtype="float32")
    if window:
        size = [max(1, int(round(window / m))) | 1 for m in (row_m, col_m)]   # odd number of cells
        height = uniform_filter(height, size, mode="nearest")
    dy, dx = np.gradient(height, row_m, col_m)
    return np.degrees(np.arctan(np.hypot(dx, dy)))


def height_above_water(height, water, step=3):
    """Height of each cell above the nearest water cell (metres); NaN if there is no water.

    A simple stand-in for "height above drainage". To save memory the nearest
    water is searched on a grid `step` times coarser, using the lowest water
    height in each coarse cell (step=3 on a 10 m grid gives 30 m, the
    resolution of the Copernicus DEM).
    """
    height = np.asarray(height, dtype="float32")
    water = np.asarray(water, dtype=bool) & np.isfinite(height)
    if not water.any():
        return np.full(height.shape, np.nan, dtype="float32")

    rows, cols = height.shape
    pr, pc = -rows % step, -cols % step
    level = np.pad(np.where(water, height, np.inf), ((0, pr), (0, pc)), constant_values=np.inf)
    level = level.reshape((rows + pr) // step, step, (cols + pc) // step, step).min(axis=(1, 3))
    has_water = np.isfinite(level)

    nearest = distance_transform_edt(~has_water, return_distances=False, return_indices=True)
    reference = level[nearest[0], nearest[1]]
    reference = np.repeat(np.repeat(reference, step, axis=0), step, axis=1)[:rows, :cols]
    return height - reference


def flood_confidence(pre_db, post_db, water_threshold=-23.3, min_drop=3.0, width=1.0, drop_width=1.0,
                     height=None, spacing=10.0, max_slope=5.0, slope_width=1.0, slope_window=90.0,
                     near_water=3.0, near_width=0.5, max_height=10.0, height_width=2.0, invalid=None):
    """Flood confidence score (0-1) and permanent-water mask.

    pre_db, post_db: images in dB on the same grid.
    width, drop_width: dB over which scores 1-3 go from 0.5 to about 0.73.
    height: optional ground height (m) on the same grid; turns on the terrain scores 5a and 5b.
        Where height is NaN, terrain is not used.
    spacing: cell size in metres, or a (row, column) pair, for the slope.
    slope_window: metres over which heights are averaged before the slope is taken.
    near_water, near_width: the slope is ignored (score 1) for cells this many metres or less above
        the nearest permanent water, with a smooth change over about near_width metres.
    invalid: optional boolean mask (e.g. layover or shadow on either date); those pixels get no score.

    Returns a dict:
        score      float32, NaN where there is no data
        permanent  bool, dark both before and after (as in change_detection)
        evidence   dict of the separate float32 scores (dark, drop, not_water_before, and
                   flat and low when height is given), for checking
    """
    pre = np.asarray(pre_db, dtype="float32")
    post = np.asarray(post_db, dtype="float32")
    if pre.shape != post.shape:
        raise ValueError(f"Pre- and post-event images differ in shape: {pre.shape} vs {post.shape}.")

    evidence = {
        "dark": rising(water_threshold - post, 0.0, width),
        "drop": rising(pre - post, min_drop, drop_width),
        "not_water_before": rising(pre - water_threshold, 0.0, width),
    }
    permanent = (post < water_threshold) & (pre < water_threshold)

    if height is not None:
        height = np.asarray(height, dtype="float32")
        if height.shape != post.shape:
            raise ValueError(f"The height grid differs in shape from the images: {height.shape} vs {post.shape}.")
        above_water = height_above_water(height, permanent)
        flat = rising(max_slope - slope_degrees(height, spacing, slope_window), 0.0, slope_width)
        flat = np.fmax(flat, rising(near_water - above_water, 0.0, near_width))   # at the water level
        low = rising(max_height - above_water, 0.0, height_width)
        del above_water
        evidence["flat"] = np.where(np.isnan(flat), np.float32(1), flat)    # no height: terrain not used
        evidence["low"] = np.where(np.isnan(low), np.float32(1), low)

    score = np.minimum.reduce(list(evidence.values()))
    nodata = np.isnan(pre) | np.isnan(post)
    if invalid is not None:
        nodata |= np.asarray(invalid, dtype=bool)
    score[nodata] = np.nan
    permanent &= ~nodata
    return {"score": score, "permanent": permanent, "evidence": evidence}


def confidence_classes(score, permanent, low=0.2, medium=0.5, high=0.8):
    """Class raster from the score: LOW from `low`, MEDIUM above `medium`, HIGH from `high`.

    MEDIUM starts just above 0.5 so that, without terrain, MEDIUM + HIGH is the
    change-detection flood. Permanent water overrides the score.
    """
    score = np.asarray(score, dtype="float64")
    out = np.full(score.shape, NOT_FLOOD, dtype="uint8")
    with np.errstate(invalid="ignore"):
        out[score >= low] = LOW
        out[score > medium] = MEDIUM
        out[score >= high] = HIGH
    out[np.asarray(permanent, dtype=bool)] = PERMANENT_WATER
    out[np.isnan(score)] = NODATA
    return out
