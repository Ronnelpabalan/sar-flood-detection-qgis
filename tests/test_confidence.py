"""Tests for core/confidence.py (no QGIS needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.change import FLOOD, PERMANENT_WATER, change_detection
from sar_flood_detection.core.confidence import (HIGH, LOW, MEDIUM, NOT_FLOOD, confidence_classes,
                                                 flood_confidence, height_above_water, rising,
                                                 slope_degrees)
from sar_flood_detection.core.thresholds import NODATA


class TestPieces(unittest.TestCase):
    def test_rising(self):
        out = rising(np.array([-2.0, 0.0, 2.0, np.nan]), 0.0, 1.0)
        np.testing.assert_allclose(out[:3], [0.1192, 0.5, 0.8808], atol=1e-4)
        self.assertTrue(np.isnan(out[3]))
        self.assertEqual(out.dtype, np.float32)

    def test_slope_of_a_ramp(self):
        height = np.tile(np.arange(6, dtype="float32") * 10.0, (4, 1))   # rises 10 m per 10 m cell
        np.testing.assert_allclose(slope_degrees(height, 10.0), 45.0, atol=1e-4)

    def test_height_above_nearest_water(self):
        height = np.array([[2.0, 3.0, 5.0, 9.0, 1.0, 1.5]], dtype="float32")
        water = np.array([[True, False, False, False, True, False]])
        np.testing.assert_allclose(height_above_water(height, water, step=1), [[0, 1, 3, 8, 0, 0.5]])

    def test_height_above_water_on_a_coarse_grid(self):
        height = np.arange(7, dtype="float32")[None, :].repeat(2, axis=0)
        water = np.zeros((2, 7), dtype=bool)
        water[:, 0] = True
        # Water only in coarse cell 0 (lowest water height 0); cells 3-6 are measured from it too.
        np.testing.assert_allclose(height_above_water(height, water, step=3), height)

    def test_no_water_gives_nan(self):
        self.assertTrue(np.isnan(height_above_water(np.zeros((3, 3)), np.zeros((3, 3), bool))).all())


class TestFloodConfidence(unittest.TestCase):
    def test_score_is_the_weakest_evidence(self):
        # post 2 dB below the threshold, drop 7 dB, pre 5 dB above the threshold: "dark" is weakest.
        r = flood_confidence(np.array([-18.3]), np.array([-25.3]))
        self.assertAlmostEqual(float(r["score"][0]), float(rising(2.0, 0.0, 1.0)), places=6)
        self.assertEqual(set(r["evidence"]), {"dark", "drop", "not_water_before"})

    def test_same_as_change_detection_without_terrain(self):
        rng = np.random.default_rng(1)
        pre = rng.uniform(-30, -5, (200, 200)).astype("float32")
        post = rng.uniform(-30, -5, (200, 200)).astype("float32")
        pre[0, 0] = np.nan
        r = flood_confidence(pre, post)
        classes = confidence_classes(r["score"], r["permanent"])
        old = change_detection(pre, post)
        np.testing.assert_array_equal(np.isin(classes, (MEDIUM, HIGH)), old == FLOOD)
        np.testing.assert_array_equal(classes == 4, old == PERMANENT_WATER)
        np.testing.assert_array_equal(classes == NODATA, old == NODATA)

    def test_terrain_only_lowers_the_score(self):
        rng = np.random.default_rng(2)
        pre = rng.uniform(-30, -5, (60, 60)).astype("float32")
        post = rng.uniform(-30, -5, (60, 60)).astype("float32")
        height = rng.uniform(0, 50, (60, 60)).astype("float32")
        plain = flood_confidence(pre, post)["score"]
        with_terrain = flood_confidence(pre, post, height=height)["score"]
        self.assertTrue(np.all(with_terrain <= plain))

    def test_hillside_is_not_flood(self):
        # Two cells that are dark after the event and dropped strongly: one at the river level on
        # flat ground, one on a steep slope high above the river.
        height = np.tile(np.array([0.0, 0.0, 0.0, 100.0, 200.0], dtype="float32"), (3, 1))
        pre = np.full((3, 5), -10.0, dtype="float32")
        post = np.full((3, 5), -28.0, dtype="float32")
        pre[:, 0] = -27.0                                       # column 0: the river (permanent water)
        r = flood_confidence(pre, post, height=height, spacing=10.0)
        classes = confidence_classes(r["score"], r["permanent"])
        self.assertTrue((classes[:, 0] == 4).all())
        self.assertTrue((classes[:, 1] == HIGH).all())        # flat, at river level
        self.assertTrue((classes[:, 4] == NOT_FLOOD).all())   # steep and 200 m up

    def test_missing_height_leaves_terrain_out(self):
        height = np.full((2, 3), np.nan, dtype="float32")
        pre, post = np.full((2, 3), -10.0), np.full((2, 3), -28.0)
        np.testing.assert_allclose(flood_confidence(pre, post, height=height)["score"],
                                   flood_confidence(pre, post)["score"])

    def test_invalid_and_nan_are_nodata(self):
        pre = np.array([[-10.0, -10.0, np.nan]])
        post = np.array([[-28.0, -28.0, -28.0]])
        r = flood_confidence(pre, post, invalid=np.array([[False, True, False]]))
        classes = confidence_classes(r["score"], r["permanent"])
        np.testing.assert_array_equal(classes, [[HIGH, NODATA, NODATA]])

    def test_shape_mismatch_raises(self):
        with self.assertRaises(ValueError):
            flood_confidence(np.zeros((2, 2)), np.zeros((2, 3)))
        with self.assertRaises(ValueError):
            flood_confidence(np.zeros((2, 2)), np.zeros((2, 2)), height=np.zeros((3, 3)))


class TestClasses(unittest.TestCase):
    def test_class_boundaries(self):
        score = np.array([0.1, 0.2, 0.5, 0.51, 0.8, np.nan, 0.9])
        permanent = np.array([False, False, False, False, False, False, True])
        np.testing.assert_array_equal(confidence_classes(score, permanent),
                                      [NOT_FLOOD, LOW, LOW, MEDIUM, HIGH, NODATA, 4])


if __name__ == "__main__":
    unittest.main()
