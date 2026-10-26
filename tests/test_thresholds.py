"""Tests for core/thresholds.py (no QGIS needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.thresholds import NODATA, classify_below, otsu_threshold, to_db


def two_groups(water_fraction, spread, rng, n=100_000, water_db=-22.0, land_db=-10.0):
    """Backscatter values from a water group and a land group with normal spread (dB)."""
    k = int(n * water_fraction)
    return np.concatenate([rng.normal(water_db, spread, k), rng.normal(land_db, spread, n - k)])


class TestToDb(unittest.TestCase):
    def test_known_values(self):
        np.testing.assert_allclose(to_db([100.0, 1.0, 0.1]), [20.0, 0.0, -10.0], atol=1e-5)

    def test_invalid_values_become_nan(self):
        self.assertTrue(np.isnan(to_db([0.0, -1.0, np.nan])).all())


class TestClassifyBelow(unittest.TestCase):
    def test_mask_values(self):
        db = np.array([[-30.0, -23.3, -10.0, np.nan]], dtype="float32")
        out = classify_below(db, -23.3)
        np.testing.assert_array_equal(out, [[1, 0, 0, NODATA]])   # -23.3 itself is not below
        self.assertEqual(out.dtype, np.uint8)


class TestOtsu(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(0)

    def test_two_clear_groups(self):
        r = otsu_threshold(two_groups(0.2, 2.0, self.rng))
        self.assertTrue(-18 < r["threshold"] < -14)
        self.assertAlmostEqual(r["below_fraction"], 0.2, delta=0.01)
        self.assertTrue(r["bimodal"])

    def test_no_water_is_flagged(self):
        r = otsu_threshold(self.rng.normal(-10.0, 2.0, 100_000))
        self.assertFalse(r["bimodal"])
        self.assertAlmostEqual(r["threshold"], -10.0, delta=0.5)   # splits the single peak

    def test_heavy_overlap_is_flagged(self):
        r = otsu_threshold(two_groups(0.2, 4.0, self.rng))
        self.assertFalse(r["bimodal"])

    def test_nan_ignored(self):
        x = two_groups(0.2, 2.0, self.rng)
        with_nan = np.concatenate([x, np.full(5000, np.nan)])
        self.assertEqual(otsu_threshold(with_nan)["threshold"], otsu_threshold(x)["threshold"])

    def test_works_on_2d_image(self):
        img = two_groups(0.25, 2.0, self.rng, n=256 * 256).reshape(256, 256)
        self.assertTrue(otsu_threshold(img)["bimodal"])

    def test_constant_image(self):
        r = otsu_threshold(np.full(1000, -12.0))
        self.assertEqual(r["threshold"], -12.0)
        self.assertFalse(r["bimodal"])

    def test_no_valid_pixels_raises(self):
        with self.assertRaises(ValueError):
            otsu_threshold(np.full(100, np.nan))


if __name__ == "__main__":
    unittest.main()
