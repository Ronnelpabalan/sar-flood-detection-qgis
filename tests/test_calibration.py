"""Tests for core/calibration.py and core/annotation.py (no QGIS needed)."""

import unittest

import numpy as np

from sar_flood_detection.core.annotation import (parse_geolocation_grid, parse_product_info,
                                                 window_for_bbox)
from sar_flood_detection.core.calibration import (NOISE_FLOOR, azimuth_noise, calibrate,
                                                  interpolate_vectors, parse_calibration, parse_noise)
from tests.s1_fixtures import annotation_xml, calibration_xml, noise_xml, simple_grid


class TestParsing(unittest.TestCase):
    def test_calibration_vectors(self):
        cal = parse_calibration(calibration_xml([(0, [0, 10], [500.0, 600.0]), (100, [0, 10], [510.0, 610.0])]))
        np.testing.assert_array_equal(cal["lines"], [0, 100])
        np.testing.assert_array_equal(cal["values"][1], [510.0, 610.0])

    def test_unknown_kind_raises(self):
        with self.assertRaises(ValueError):
            parse_calibration(calibration_xml([(0, [0], [1.0])]), "sigma0")

    def test_noise_new_format(self):
        n = parse_noise(noise_xml([(0, [0, 10], [5.0, 6.0])], [(0, 99, 0, 4, [0, 99], [1.0, 2.0])]))
        self.assertEqual(len(n["azimuth"]), 1)
        self.assertEqual(n["azimuth"][0]["last_sample"], 4)

    def test_noise_old_format_has_no_azimuth_table(self):
        n = parse_noise(noise_xml([(0, [0, 10], [5.0, 6.0])], old_format=True))
        self.assertEqual(n["azimuth"], [])
        np.testing.assert_array_equal(n["range"]["values"][0], [5.0, 6.0])


class TestInterpolation(unittest.TestCase):
    def setUp(self):
        # value = 1000 + line + 10 * pixel, given every 100 lines and 50 pixels
        px = [0, 50, 100]
        self.vectors = parse_calibration(calibration_xml(
            [(ln, px, [1000 + ln + 10 * p for p in px]) for ln in (0, 100, 200)]))

    def test_exact_at_vector_points(self):
        out = interpolate_vectors(self.vectors, [0, 100, 200], [0, 50, 100])
        expected = [[1000 + ln + 10 * p for p in (0, 50, 100)] for ln in (0, 100, 200)]
        np.testing.assert_allclose(out, expected)

    def test_bilinear_between_points(self):
        rows, cols = np.arange(30, 170, 7), np.arange(3, 97, 11)
        out = interpolate_vectors(self.vectors, rows, cols)
        np.testing.assert_allclose(out, 1000 + rows[:, None] + 10 * cols[None, :])

    def test_nearest_value_outside_vectors(self):
        out = interpolate_vectors(self.vectors, [250], [120])
        self.assertAlmostEqual(out[0, 0], 1000 + 200 + 10 * 100)

    def test_single_vector(self):
        one = parse_calibration(calibration_xml([(0, [0, 10], [1.0, 2.0])]))
        np.testing.assert_allclose(interpolate_vectors(one, [5, 6], [0, 5, 10]), [[1, 1.5, 2]] * 2)

    def test_vectors_with_different_pixel_positions(self):
        v = parse_noise(noise_xml([(0, [0, 100], [0.0, 100.0]), (10, [0, 40, 100], [0.0, 40.0, 100.0])]))
        out = interpolate_vectors(v["range"], [0, 5, 10], [20, 70])
        np.testing.assert_allclose(out, [[20, 70]] * 3)

    def test_azimuth_blocks_apply_to_their_samples_only(self):
        blocks = parse_noise(noise_xml([(0, [0], [1.0])], [(0, 9, 0, 4, [0, 9], [2.0, 2.0]),
                                                            (0, 9, 5, 9, [0, 9], [1.0, 10.0])]))["azimuth"]
        out = azimuth_noise(blocks, np.arange(10), np.arange(12))
        np.testing.assert_allclose(out[:, :5], 2.0)
        np.testing.assert_allclose(out[:, 5:10], np.arange(1.0, 11.0)[:, None] * np.ones(5))
        np.testing.assert_allclose(out[:, 10:], 1.0)         # beyond the last block


class TestCalibrate(unittest.TestCase):
    def setUp(self):
        self.cal = parse_calibration(calibration_xml([(0, [0, 100], [500.0, 500.0]),
                                                      (100, [0, 100], [500.0, 500.0])]))

    def test_formula_without_noise(self):
        dn = np.full((4, 5), 50.0)
        np.testing.assert_allclose(calibrate(dn, 0, 0, self.cal), (50.0 / 500.0) ** 2, rtol=1e-6)

    def test_noise_is_subtracted_before_calibration(self):
        noise = parse_noise(noise_xml([(0, [0, 100], [900.0, 900.0])],
                                      [(0, 100, 0, 100, [0, 100], [1.0, 1.0])]))
        dn = np.full((3, 3), 50.0)                            # 2500 - 900 = 1600 = 40**2
        np.testing.assert_allclose(calibrate(dn, 10, 20, self.cal, noise), (40.0 / 500.0) ** 2, rtol=1e-6)

    def test_azimuth_factor_scales_range_noise(self):
        noise = parse_noise(noise_xml([(0, [0, 100], [450.0, 450.0])],
                                      [(0, 100, 0, 100, [0, 100], [2.0, 2.0])]))
        np.testing.assert_allclose(calibrate(np.full((2, 2), 50.0), 0, 0, self.cal, noise),
                                   (40.0 / 500.0) ** 2, rtol=1e-6)

    def test_negative_after_noise_removal_gets_floor(self):
        noise = parse_noise(noise_xml([(0, [0, 100], [5000.0, 5000.0])], old_format=True))
        out = calibrate(np.full((2, 2), 50.0), 0, 0, self.cal, noise)
        np.testing.assert_allclose(out, NOISE_FLOOR)

    def test_zero_dn_is_no_data(self):
        dn = np.array([[0.0, 50.0]])
        out = calibrate(dn, 0, 0, self.cal)
        self.assertTrue(np.isnan(out[0, 0]))
        self.assertFalse(np.isnan(out[0, 1]))
        self.assertEqual(out.dtype, np.float32)

    def test_window_offset_matches_full_image(self):
        cal = parse_calibration(calibration_xml([(ln, [0, 60], [400.0 + ln, 700.0 + ln]) for ln in (0, 40, 80)]))
        noise = parse_noise(noise_xml([(ln, [0, 60], [100.0 + ln, 300.0]) for ln in (0, 40, 80)],
                                      [(0, 80, 0, 29, [0, 80], [1.0, 1.5]), (0, 80, 30, 60, [0, 80], [0.8, 1.2])]))
        dn = np.random.default_rng(0).uniform(20, 200, (80, 60))
        full = calibrate(dn, 0, 0, cal, noise)
        part = calibrate(dn[25:55, 10:45], 25, 10, cal, noise)
        np.testing.assert_allclose(part, full[25:55, 10:45], rtol=1e-6)


class TestAnnotation(unittest.TestCase):
    def setUp(self):
        self.xml = annotation_xml(400, 800, simple_grid(400, 800))

    def test_product_info(self):
        info = parse_product_info(self.xml)
        self.assertEqual((info["mission"], info["polarisation"], info["pass"]), ("S1D", "VH", "Descending"))
        self.assertEqual((info["lines"], info["samples"]), (400, 800))

    def test_geolocation_grid(self):
        g = parse_geolocation_grid(self.xml)
        self.assertEqual(len(g["line"]), 25)
        self.assertAlmostEqual(g["longitude"][-1], 105.8)

    def test_window_for_box_inside_image(self):
        g = parse_geolocation_grid(self.xml)
        # lon 105.1..105.3 -> pixels 100..300; lat 19.9..19.85 -> lines 100..150
        self.assertEqual(window_for_bbox(g, (105.1, 19.85, 105.3, 19.9), 400, 800, margin=10),
                         (90, 90, 220, 70))

    def test_window_is_clipped_to_image(self):
        g = parse_geolocation_grid(self.xml)
        self.assertEqual(window_for_bbox(g, (104.0, 19.0, 106.0, 21.0), 400, 800, margin=10), (0, 0, 800, 400))

    def test_box_outside_image_raises(self):
        g = parse_geolocation_grid(self.xml)
        with self.assertRaises(ValueError):
            window_for_bbox(g, (110.0, 10.0, 111.0, 11.0), 400, 800)


if __name__ == "__main__":
    unittest.main()
