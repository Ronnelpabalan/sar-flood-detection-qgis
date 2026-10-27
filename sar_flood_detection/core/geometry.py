"""Range-Doppler geometry: where a point on the ground appears in a Sentinel-1 GRD image.

For a ground point the satellite sees it at the moment its line of sight is
perpendicular to the flight direction (zero Doppler). That time gives the
image line; the distance at that time (slant range) gives the image column,
after ESA's slant-to-ground-range polynomials. This is the same model that
SNAP's Range-Doppler Terrain Correction uses.
"""

import xml.etree.ElementTree as ET

import numpy as np

from .orbit import parse_time

SPEED_OF_LIGHT = 299792458.0                       # m/s
WGS84_A = 6378137.0                                # semi-major axis (m)
WGS84_E2 = 6.69437999014e-3                        # first eccentricity squared


def geodetic_to_ecef(lat, lon, height):
    """WGS84 latitude, longitude (degrees) and ellipsoidal height (m) -> Earth-fixed x, y, z (m), shape (..., 3)."""
    lat = np.radians(np.asarray(lat, dtype="float64"))
    lon = np.radians(np.asarray(lon, dtype="float64"))
    height = np.asarray(height, dtype="float64")
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    return np.stack([(n + height) * np.cos(lat) * np.cos(lon),
                     (n + height) * np.cos(lat) * np.sin(lon),
                     (n * (1.0 - WGS84_E2) + height) * np.sin(lat)], -1)


def parse_geometry(xml_text):
    """Image timing and slant-to-ground-range polynomials from annotation/s1*.xml."""
    root = ET.fromstring(xml_text)
    image = root.find("imageAnnotation/imageInformation")
    times, sr0, coefficients = [], [], []
    for c in root.findall(".//coordinateConversionList/coordinateConversion"):   # the parent has the same tag
        times.append(parse_time(c.findtext("azimuthTime")))
        sr0.append(float(c.findtext("sr0")))
        coefficients.append([float(v) for v in c.findtext("srgrCoefficients").split()])
    if not times:
        raise ValueError("No coordinateConversion elements found (not a GRD product?).")
    order = np.argsort(times)
    return {
        "first_line_time": parse_time(image.findtext("productFirstLineUtcTime")),
        "line_interval": float(image.findtext("azimuthTimeInterval")),
        "near_range": float(image.findtext("slantRangeTime")) * SPEED_OF_LIGHT / 2,
        "range_spacing": float(image.findtext("rangePixelSpacing")),
        "lines": int(image.findtext("numberOfLines")),
        "samples": int(image.findtext("numberOfSamples")),
        "wavelength": SPEED_OF_LIGHT / float(root.findtext(".//productInformation/radarFrequency")),
        "srgr_times": np.array(times)[order],
        "sr0": np.array(sr0)[order],
        "srgr": np.array(coefficients)[order],
    }


def _polyval(coefficients, x):
    """Per-point polynomial c0 + c1 x + c2 x² + ... (coefficients shape (n, k), x shape (n,))."""
    result = np.zeros_like(x)
    for i in range(coefficients.shape[1] - 1, -1, -1):
        result = result * x + coefficients[:, i]
    return result


class RadarGeometry:
    """Maps ground points to (line, pixel) in one GRD image, given its orbit."""

    def __init__(self, params, orbit):
        self.p, self.orbit = params, orbit

    def zero_doppler(self, xyz, iterations=20, tolerance=1e-8):
        """Zero-Doppler time (s since 1970) and slant range (m) for Earth-fixed points xyz (n, 3).

        Newton's method on f(t) = (P(t) - X) . V(t), whose derivative is
        V . V + (P - X) . A. It starts at the middle of the image.
        """
        xyz = np.atleast_2d(np.asarray(xyz, dtype="float64"))
        t_min, t_max = self.orbit.times[0], self.orbit.times[-1]
        mid = self.p["first_line_time"] + self.p["line_interval"] * (self.p["lines"] - 1) / 2
        t = np.full(len(xyz), mid)
        for _ in range(iterations):
            d = self.orbit.position(t) - xyz
            v = self.orbit.velocity(t)
            f = np.sum(d * v, 1)
            df = np.sum(v * v, 1) + np.sum(d * self.orbit.acceleration(t), 1)
            step = f / df
            t = np.clip(t - step, t_min, t_max)
            if np.abs(step).max() < tolerance:
                break
        return t, np.linalg.norm(self.orbit.position(t) - xyz, axis=1)

    def ground_range(self, slant_range, t):
        """ESA's slant-to-ground-range polynomials, linearly interpolated in time between records."""
        times = self.p["srgr_times"]
        if len(times) == 1:
            return _polyval(self.p["srgr"][[0] * len(t)], slant_range - self.p["sr0"][0])
        i = np.clip(np.searchsorted(times, t), 1, len(times) - 1)
        w = np.clip((t - times[i - 1]) / (times[i] - times[i - 1]), 0.0, 1.0)
        g0 = _polyval(self.p["srgr"][i - 1], slant_range - self.p["sr0"][i - 1])
        g1 = _polyval(self.p["srgr"][i], slant_range - self.p["sr0"][i])
        return (1 - w) * g0 + w * g1

    def locate(self, xyz):
        """(line, pixel) as float arrays for Earth-fixed points xyz (n, 3); may fall outside the image."""
        t, slant_range = self.zero_doppler(xyz)
        line = (t - self.p["first_line_time"]) / self.p["line_interval"]
        pixel = self.ground_range(slant_range, t) / self.p["range_spacing"]
        return line, pixel

    def to_image(self, lat, lon, height):
        """(line, pixel) for WGS84 latitude, longitude (degrees) and ellipsoidal height (m)."""
        xyz = geodetic_to_ecef(np.ravel(lat), np.ravel(lon), np.ravel(height))
        line, pixel = self.locate(xyz)
        shape = np.shape(lat)
        return line.reshape(shape), pixel.reshape(shape)
