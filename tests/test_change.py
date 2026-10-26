"""Tests for core/change.py (no QGIS needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.change import DRY, FLOOD, PERMANENT_WATER, change_detection
from sar_flood_detection.core.thresholds import NODATA


class TestChangeDetection(unittest.TestCase):
    def test_each_class(self):
        #                 land->flood  river  land  dark land, small drop  missing
        pre = np.array([[-10.0,        -26.0, -10.0, -22.0,                  np.nan]])
        post = np.array([[-25.0,       -27.0, -11.0, -24.0,                  -25.0]])
        out = change_detection(pre, post, water_threshold=-23.3, min_drop=3.0)
        np.testing.assert_array_equal(out, [[FLOOD, PERMANENT_WATER, DRY, DRY, NODATA]])
        self.assertEqual(out.dtype, np.uint8)

    def test_drop_must_exceed_min_drop(self):
        pre, post = np.array([-21.0]), np.array([-24.0])            # drop of exactly 3 dB
        self.assertEqual(change_detection(pre, post, min_drop=3.0)[0], DRY)
        self.assertEqual(change_detection(pre, post, min_drop=2.9)[0], FLOOD)

    def test_nan_after_event_is_nodata(self):
        out = change_detection(np.array([-10.0]), np.array([np.nan]))
        self.assertEqual(out[0], NODATA)

    def test_shape_mismatch_raises(self):
        with self.assertRaises(ValueError):
            change_detection(np.zeros((2, 2)), np.zeros((2, 3)))


if __name__ == "__main__":
    unittest.main()
