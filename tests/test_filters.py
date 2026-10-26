"""Tests for core/filters.py on synthetic speckled images (no QGIS needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.filters import lee_filter_db, median_filter_db

ENL = 4.4


def speckle(true_db, rng):
    """Multiply the true backscatter by gamma-distributed speckle with mean 1."""
    linear = 10 ** (true_db / 10) * rng.gamma(ENL, 1 / ENL, true_db.shape)
    return (10 * np.log10(linear)).astype("float32")


class TestFilters(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(0)
        self.filters = {"lee": lee_filter_db, "median": median_filter_db}

    def test_constant_image_unchanged(self):
        img = np.full((64, 64), -12.0, "float32")
        for name, f in self.filters.items():
            with self.subTest(name):
                np.testing.assert_allclose(f(img), -12.0, atol=1e-4)

    def test_reduces_speckle_on_flat_area(self):
        img = speckle(np.full((256, 256), -10.0), self.rng)
        for name, f in self.filters.items():
            with self.subTest(name):
                self.assertLess(f(img).std(), 0.6 * img.std())

    def test_lee_preserves_mean_power(self):
        img = speckle(np.full((256, 256), -10.0), self.rng)
        out = lee_filter_db(img)
        self.assertAlmostEqual(np.mean(10 ** (out / 10)), np.mean(10 ** (img / 10)), places=3)

    def test_fewer_errors_at_water_land_edge(self):
        truth = np.zeros((256, 256), bool)
        truth[:, :128] = True
        img = speckle(np.where(truth, -20.0, -8.0), self.rng)
        raw_errors = np.mean((img < -14) != truth)
        for name, f in self.filters.items():
            with self.subTest(name):
                self.assertLess(np.mean((f(img) < -14) != truth), raw_errors)

    def test_nan_pixels_stay_nan_and_neighbours_finite(self):
        img = speckle(np.full((64, 64), -10.0), self.rng)
        img[20:25, 20:25] = np.nan
        for name, f in self.filters.items():
            with self.subTest(name):
                out = f(img)
                self.assertTrue(np.isnan(out[20:25, 20:25]).all())
                self.assertEqual(np.isfinite(out).sum(), 64 * 64 - 25)

    def test_all_nan_image(self):
        img = np.full((16, 16), np.nan, "float32")
        for name, f in self.filters.items():
            with self.subTest(name):
                self.assertTrue(np.isnan(f(img)).all())

    def test_larger_window_smooths_more(self):
        img = speckle(np.full((256, 256), -10.0), self.rng)
        for name, f in self.filters.items():
            with self.subTest(name):
                self.assertLess(f(img, size=9).std(), f(img, size=3).std())


if __name__ == "__main__":
    unittest.main()
