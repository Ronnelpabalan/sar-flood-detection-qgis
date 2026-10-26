"""Thresholds that separate water from land in SAR backscatter (dB).

Masks are uint8: 1 = water, 0 = not water, NODATA = no valid input.
"""

import numpy as np

NODATA = 255


def to_db(linear):
    """Convert linear power to dB. Zero, negative and NaN values become NaN."""
    linear = np.asarray(linear, dtype="float64")
    with np.errstate(divide="ignore", invalid="ignore"):
        db = 10.0 * np.log10(np.where(linear > 0, linear, np.nan))
    return db.astype("float32")


def classify_below(db, threshold):
    """Water mask: 1 where db < threshold, 0 elsewhere, NODATA where db is NaN."""
    out = np.where(db < threshold, 1, 0).astype("uint8")
    out[np.isnan(db)] = NODATA
    return out


def otsu_threshold(db, nbins=256, max_valley_ratio=0.5):
    """Otsu threshold of the valid values, with a check that the histogram is two-peaked.

    Otsu picks the value that best splits the histogram into two groups. It
    assumes two peaks (water and land); with little water in the image it
    splits the land peak instead. The check smooths the histogram, finds the
    highest point on each side of the threshold, and compares the lowest point
    between them (the valley) with the lower of the two peaks. A deep valley
    means two real peaks; a ratio near 1 means one peak was cut in half.

    Returns a dict:
        threshold       Otsu threshold (dB)
        below_fraction  share of valid pixels below the threshold
        valley_ratio    valley height / lower peak height (0 = clear valley, 1 = none)
        bimodal         valley_ratio <= max_valley_ratio (a heuristic, not a test)
    """
    x = np.asarray(db, dtype="float64")
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("The image has no valid pixels.")
    lo, hi = np.percentile(x, [0.1, 99.9])          # ignore extreme speckle values when binning
    if hi <= lo:
        return {"threshold": float(lo), "below_fraction": 0.0, "valley_ratio": 1.0, "bimodal": False}

    counts, edges = np.histogram(np.clip(x, lo, hi), bins=nbins, range=(lo, hi))
    centres = (edges[:-1] + edges[1:]) / 2
    p = counts / counts.sum()
    w0 = np.cumsum(p)                                # weight of the class below each split
    mu = np.cumsum(p * centres)
    w1 = 1.0 - w0
    with np.errstate(divide="ignore", invalid="ignore"):
        between = (mu[-1] * w0 - mu) ** 2 / (w0 * w1)   # between-class variance
    between[~np.isfinite(between)] = 0.0
    k = min(int(np.argmax(between)), nbins - 2)      # split between bins k and k + 1
    threshold = float(edges[k + 1])

    smooth = np.convolve(counts, np.ones(5) / 5, mode="same")
    i_low = int(np.argmax(smooth[: k + 1]))
    i_high = k + 1 + int(np.argmax(smooth[k + 1:]))
    lower_peak = min(smooth[i_low], smooth[i_high])
    valley = smooth[i_low: i_high + 1].min()
    valley_ratio = float(valley / lower_peak) if lower_peak > 0 else 1.0

    return {
        "threshold": threshold,
        "below_fraction": float(np.mean(x < threshold)),
        "valley_ratio": valley_ratio,
        "bimodal": valley_ratio <= max_valley_ratio,
    }
