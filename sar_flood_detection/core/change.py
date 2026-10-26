"""Flood mapping by change detection between a pre-event and a post-event image.

Both images must be in dB and on the same pixel grid.
Classes: DRY = 0, FLOOD = 1, PERMANENT_WATER = 2, NODATA = 255.
"""

import numpy as np

from .thresholds import NODATA

DRY = 0
FLOOD = 1
PERMANENT_WATER = 2


def change_detection(pre_db, post_db, water_threshold=-23.3, min_drop=3.0):
    """Classify each pixel by comparing the pre- and post-event backscatter.

    PERMANENT_WATER: below water_threshold both before and after the event.
    FLOOD:           below water_threshold after the event, not before, and
                     backscatter dropped by more than min_drop dB.
    DRY:             everything else, including pixels that are dark after the
                     event but changed little (uncertain, kept conservative).
    NODATA:          either image has no value.
    """
    pre = np.asarray(pre_db, dtype="float32")
    post = np.asarray(post_db, dtype="float32")
    if pre.shape != post.shape:
        raise ValueError(f"Pre- and post-event images differ in shape: {pre.shape} vs {post.shape}.")

    water_after = post < water_threshold
    water_before = pre < water_threshold
    dropped = (pre - post) > min_drop

    out = np.full(post.shape, DRY, dtype="uint8")
    out[water_after & water_before] = PERMANENT_WATER
    out[water_after & ~water_before & dropped] = FLOOD
    out[np.isnan(pre) | np.isnan(post)] = NODATA
    return out
