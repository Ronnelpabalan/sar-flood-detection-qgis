"""Radiometric calibration and thermal noise removal for Sentinel-1 GRD.

Follows ESA's Sentinel-1 calibration and noise removal documents:

    sigma0 = (DN**2 - noise) / A**2

DN is the GRD pixel value, A the calibration value (sigmaNought, betaNought
or gamma) and noise the thermal noise power, both read from the XML files in
the product's annotation/calibration folder. The look-up tables are given at
sparse lines and pixels and are interpolated bilinearly to every pixel.

Since IPF 2.9 (2018) the noise is split into a range table and an azimuth
table that are multiplied; older products have only a range table.

Arrays are handled in windows: rows and cols are the image line and pixel
numbers of the window, so a small area can be calibrated without reading
the whole scene.
"""

import xml.etree.ElementTree as ET

import numpy as np

CALIBRATION_KINDS = ("sigmaNought", "betaNought", "gamma")
NOISE_FLOOR = 1e-5          # linear power (-50 dB) used where noise removal gives zero or less


def _numbers(text):
    return np.array(text.split(), dtype="float64")


def parse_calibration(xml_text, kind="sigmaNought"):
    """Calibration vectors from a calibration-*.xml file.

    Returns {"lines": array, "pixels": [array per vector], "values": [array per vector]}.
    """
    if kind not in CALIBRATION_KINDS:
        raise ValueError(f"kind must be one of {CALIBRATION_KINDS}, got {kind!r}")
    root = ET.fromstring(xml_text)
    lines, pixels, values = [], [], []
    for v in root.iter("calibrationVector"):
        lines.append(int(v.findtext("line")))
        pixels.append(_numbers(v.findtext("pixel")))
        values.append(_numbers(v.findtext(kind)))
    if not lines:
        raise ValueError("No calibrationVector elements found.")
    return {"lines": np.array(lines), "pixels": pixels, "values": values}


def parse_noise(xml_text):
    """Thermal noise tables from a noise-*.xml file.

    Returns {"range": {"lines", "pixels", "values"}, "azimuth": [block, ...]},
    where each azimuth block is a dict with first_line, last_line,
    first_sample, last_sample, lines and values. "azimuth" is empty for
    products made before IPF 2.9.
    """
    root = ET.fromstring(xml_text)
    lines, pixels, values = [], [], []
    vectors = list(root.iter("noiseRangeVector"))
    lut_tag = "noiseRangeLut"
    if not vectors:                                  # products before IPF 2.9
        vectors, lut_tag = list(root.iter("noiseVector")), "noiseLut"
    for v in vectors:
        lines.append(int(v.findtext("line")))
        pixels.append(_numbers(v.findtext("pixel")))
        values.append(_numbers(v.findtext(lut_tag)))
    if not lines:
        raise ValueError("No noise range vectors found.")

    azimuth = []
    for v in root.iter("noiseAzimuthVector"):
        azimuth.append({
            "first_line": int(v.findtext("firstAzimuthLine")),
            "last_line": int(v.findtext("lastAzimuthLine")),
            "first_sample": int(v.findtext("firstRangeSample")),
            "last_sample": int(v.findtext("lastRangeSample")),
            "lines": _numbers(v.findtext("line")),
            "values": _numbers(v.findtext("noiseAzimuthLut")),
        })
    return {"range": {"lines": np.array(lines), "pixels": pixels, "values": values}, "azimuth": azimuth}


def interpolate_vectors(vectors, rows, cols):
    """Bilinear interpolation of vectors given at sparse lines and pixels.

    Each vector is interpolated along the pixels to the window columns, then
    the result is interpolated between vector lines. Outside the first and
    last vector the nearest value is used. Returns a (len(rows), len(cols))
    float64 array.
    """
    rows = np.asarray(rows, dtype="float64")
    cols = np.asarray(cols, dtype="float64")
    lines = vectors["lines"]
    # only the vectors around the window are needed
    first = max(int(np.searchsorted(lines, rows.min(), side="right")) - 1, 0)
    last = min(int(np.searchsorted(lines, rows.max(), side="left")), len(lines) - 1)
    idx = range(first, last + 1)
    at_cols = np.stack([np.interp(cols, vectors["pixels"][i], vectors["values"][i]) for i in idx])
    used_lines = lines[first:last + 1].astype("float64")
    if len(used_lines) == 1:
        return np.repeat(at_cols, len(rows), axis=0)

    k = np.clip(np.searchsorted(used_lines, rows, side="right") - 1, 0, len(used_lines) - 2)
    t = np.clip((rows - used_lines[k]) / (used_lines[k + 1] - used_lines[k]), 0.0, 1.0)
    return at_cols[k] * (1.0 - t)[:, None] + at_cols[k + 1] * t[:, None]


def azimuth_noise(blocks, rows, cols):
    """Azimuth noise factor for the window: 1 where no block applies."""
    rows = np.asarray(rows)
    cols = np.asarray(cols)
    out = np.ones((len(rows), len(cols)))
    for b in blocks:
        r = (rows >= b["first_line"]) & (rows <= b["last_line"])
        c = (cols >= b["first_sample"]) & (cols <= b["last_sample"])
        if r.any() and c.any():
            factor = np.interp(rows[r], b["lines"], b["values"])
            out[np.ix_(r, c)] = factor[:, None]
    return out


def calibrate(dn, row0, col0, calibration, noise=None, floor=NOISE_FLOOR):
    """Calibrated backscatter (linear power, float32) for a window of GRD pixel values.

    dn          2-D array of GRD values (DN) for the window
    row0, col0  image line and pixel of the window's top-left corner
    calibration output of parse_calibration (its kind sets sigma0, beta0 or gamma0)
    noise       output of parse_noise, or None to skip thermal noise removal

    DN = 0 marks no data (e.g. removed border noise) and becomes NaN. Where
    noise removal leaves zero or less, the value is set to floor.
    """
    dn = np.asarray(dn, dtype="float64")
    rows = np.arange(row0, row0 + dn.shape[0])
    cols = np.arange(col0, col0 + dn.shape[1])
    power = dn ** 2
    if noise is not None:
        noise_power = interpolate_vectors(noise["range"], rows, cols)
        if noise["azimuth"]:
            noise_power *= azimuth_noise(noise["azimuth"], rows, cols)
        power = power - noise_power
    a = interpolate_vectors(calibration, rows, cols)
    out = power / a ** 2
    out = np.where(out > 0, out, floor)
    out[dn == 0] = np.nan
    return out.astype("float32")
