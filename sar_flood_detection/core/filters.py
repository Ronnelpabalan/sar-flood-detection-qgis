"""Speckle filters for Sentinel-1 backscatter in dB, tolerant of NaN pixels.

Both filters take and return a 2-D array in dB. NaN pixels stay NaN and are
not used when filtering their neighbours (Lee) or are temporarily filled (median).
Adapted from github.com/Ronnelpabalan/sar-flood-unet (src/filters.py).
"""

import numpy as np
from scipy.ndimage import median_filter, uniform_filter


def _local_mean(x, valid, size):
    """Mean of x over the valid pixels in a size x size window (NaN if there are none)."""
    num = uniform_filter(np.where(valid, x, 0.0), size, mode="reflect")
    den = uniform_filter(valid.astype("float64"), size, mode="reflect")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 1e-6, num / den, np.nan)


def lee_filter_db(db, size=5, enl=4.4):
    """Lee speckle filter, applied in linear power; input and output in dB.

    Each pixel is pulled towards its local mean, more strongly in flat areas
    (where local variation is close to what speckle alone would cause) and
    less near edges.

    size: window size in pixels
    enl:  equivalent number of looks of the data; about 4.4 for Sentinel-1 IW GRDH
    """
    valid = ~np.isnan(db)
    x = np.where(valid, 10.0 ** (db / 10.0), 0.0)
    mean = _local_mean(x, valid, size)
    var = np.maximum(_local_mean(x ** 2, valid, size) - mean ** 2, 0.0)
    cu2 = 1.0 / enl                                            # speckle variance (relative)
    var_signal = np.maximum((var - cu2 * mean ** 2) / (1.0 + cu2), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.clip(np.where(var > 0, var_signal / var, 0.0), 0.0, 1.0)
        out = 10.0 * np.log10(np.maximum(mean + w * (x - mean), 1e-10))
    out[~valid] = np.nan
    return out.astype("float32")


def median_filter_db(db, size=5):
    """Median filter in dB (the median is the same in dB and linear units).

    NaN pixels are temporarily filled with the chip median so they do not
    disturb the filter, then set back to NaN.
    """
    valid = ~np.isnan(db)
    if not valid.any():
        return db.astype("float32", copy=True)
    filled = np.where(valid, db, np.nanmedian(db))
    out = median_filter(filled, size=size, mode="reflect")
    out[~valid] = np.nan
    return out.astype("float32")
