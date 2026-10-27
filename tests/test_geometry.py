"""Tests for core/geometry.py on a synthetic orbit (no QGIS or real data needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.geometry import (SPEED_OF_LIGHT, WGS84_A, RadarGeometry, geodetic_to_ecef,
                                               parse_geometry)
from sar_flood_detection.core.orbit import Orbit, parse_eof
from tests.s1_fixtures import eof_xml, geometry_xml
from tests.test_orbit import T0, circular, utc

FIRST, INTERVAL, LINES = 100.0, 0.01, 2500           # image: 100 s to 125 s after T0
SR0 = (799000.0, 799100.0)                           # near range of the two srgr records
SRGR_TIMES = (95.0, 130.0)


def make_geometry():
    steps = np.arange(0, 301, 10.0)
    p, v = circular(steps)
    orbit = Orbit(parse_eof(eof_xml([utc(s) for s in steps], p, v)), T0.timestamp() + FIRST,
                  T0.timestamp() + FIRST + INTERVAL * LINES)
    xml = geometry_xml(utc(FIRST), INTERVAL, 2 * SR0[0] / SPEED_OF_LIGHT, 10.0, LINES, 20000, 5.405e9,
                       [(utc(SRGR_TIMES[0]), SR0[0], [0.0, 1.0]),       # ground range = slant range - sr0
                        (utc(SRGR_TIMES[1]), SR0[1], [0.0, 1.0])])
    return RadarGeometry(parse_geometry(xml), orbit), orbit


def targets(orbit, seconds, ranges, look_angle=np.radians(35.0)):
    """Points seen at exactly these times (s after T0) and slant ranges: perpendicular to the velocity."""
    t = T0.timestamp() + np.asarray(seconds)
    p, v = orbit.position(t), orbit.velocity(t)
    along = v / np.linalg.norm(v, axis=1, keepdims=True)
    up = p - np.sum(p * along, 1, keepdims=True) * along
    down = -up / np.linalg.norm(up, axis=1, keepdims=True)
    side = np.cross(along, down)
    look = np.cos(look_angle) * down + np.sin(look_angle) * side
    return p + np.asarray(ranges)[:, None] * look


class TestGeodetic(unittest.TestCase):
    def test_known_points(self):
        np.testing.assert_allclose(geodetic_to_ecef(0, 0, 0), [WGS84_A, 0, 0], atol=1e-6)
        np.testing.assert_allclose(geodetic_to_ecef(0, 90, 100), [0, WGS84_A + 100, 0], atol=1e-6)
        np.testing.assert_allclose(geodetic_to_ecef(90, 0, 0), [0, 0, 6356752.314], atol=1e-3)   # polar radius


class TestParseGeometry(unittest.TestCase):
    def test_fields(self):
        g, _ = make_geometry()
        self.assertAlmostEqual(g.p["first_line_time"] - T0.timestamp(), FIRST, places=5)
        self.assertEqual((g.p["line_interval"], g.p["lines"], g.p["range_spacing"]), (INTERVAL, LINES, 10.0))
        self.assertAlmostEqual(g.p["near_range"], SR0[0], places=3)
        self.assertAlmostEqual(g.p["wavelength"], 0.05547, places=5)
        np.testing.assert_allclose(g.p["sr0"], SR0)                  # only the records, not the parent element
        self.assertEqual(g.p["srgr"].shape, (2, 2))


class TestRangeDoppler(unittest.TestCase):
    def setUp(self):
        self.geo, self.orbit = make_geometry()
        rng = np.random.default_rng(1)
        self.seconds = rng.uniform(FIRST, FIRST + INTERVAL * LINES, 50)
        self.ranges = rng.uniform(800e3, 950e3, 50)
        self.xyz = targets(self.orbit, self.seconds, self.ranges)

    def test_zero_doppler_time_and_range(self):
        t, r = self.geo.zero_doppler(self.xyz)
        self.assertLess(np.abs(t - T0.timestamp() - self.seconds).max(), 1e-6)     # seconds
        self.assertLess(np.abs(r - self.ranges).max(), 1e-3)                       # metres

    def test_line_and_pixel(self):
        line, pixel = self.geo.locate(self.xyz)
        w = (self.seconds - SRGR_TIMES[0]) / (SRGR_TIMES[1] - SRGR_TIMES[0])
        sr0 = (1 - w) * SR0[0] + w * SR0[1]                          # srgr records interpolated in time
        np.testing.assert_allclose(line, (self.seconds - FIRST) / INTERVAL, atol=1e-3)
        np.testing.assert_allclose(pixel, (self.ranges - sr0) / 10.0, atol=1e-3)

    def test_to_image_keeps_shape(self):
        lat, lon = np.meshgrid(np.linspace(19.6, 20.2, 3), np.linspace(105.4, 106.0, 4))
        line, pixel = self.geo.to_image(lat, lon, np.zeros_like(lat))
        self.assertEqual(line.shape, (4, 3))
        xyz = geodetic_to_ecef(lat.ravel(), lon.ravel(), 0.0 * lat.ravel())
        np.testing.assert_allclose(pixel.ravel(), self.geo.locate(xyz)[1])


if __name__ == "__main__":
    unittest.main()
