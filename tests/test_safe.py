"""Tests for safe.py on a tiny fake SAFE product (needs GDAL, e.g. QGIS's Python)."""

import os
import tempfile
import unittest
import zipfile

import numpy as np

try:
    from osgeo import gdal
    HAVE_GDAL = True
except ImportError:
    HAVE_GDAL = False

from tests.s1_fixtures import annotation_xml, calibration_xml, noise_xml, simple_grid

LINES, SAMPLES = 60, 80
SAFE = "S1D_IW_GRDH_1SDV_20261003T225043_20261003T225108_004857_009216_0529.SAFE"


def build_safe(folder, as_zip=True):
    """Write a two-polarisation fake product; returns its path and the DN arrays."""
    rng = np.random.default_rng(0)
    dn = {"VH": rng.integers(1, 300, (LINES, SAMPLES)).astype("uint16"),
          "VV": rng.integers(1, 900, (LINES, SAMPLES)).astype("uint16")}
    dn["VH"][:, :3] = 0                                         # border noise already removed
    files = {f"{SAFE}/manifest.safe": "<manifest/>"}
    tiffs = {}
    for i, pol in enumerate(("VV", "VH"), start=1):
        stem = f"s1d-iw-grd-{pol.lower()}-20261003t225043-20261003t225108-004857-009216-00{i}"
        files[f"{SAFE}/annotation/{stem}.xml"] = annotation_xml(LINES, SAMPLES, simple_grid(LINES, SAMPLES), pol)
        files[f"{SAFE}/annotation/calibration/calibration-{stem}.xml"] = calibration_xml(
            [(0, [0, SAMPLES], [500.0, 600.0]), (LINES, [0, SAMPLES], [520.0, 620.0])])
        files[f"{SAFE}/annotation/calibration/noise-{stem}.xml"] = noise_xml(
            [(0, [0, SAMPLES], [100.0, 50.0]), (LINES, [0, SAMPLES], [100.0, 50.0])],
            [(0, LINES - 1, 0, SAMPLES - 1, [0, LINES - 1], [1.0, 1.0])])
        tiffs[f"{SAFE}/measurement/{stem}.tiff"] = dn[pol]
    files[f"{SAFE}/support/s1-level-1-calibration.xsd"] = "<xsd/>"     # must be ignored

    staging = os.path.join(folder, "staging")
    for name, content in files.items():
        os.makedirs(os.path.dirname(os.path.join(staging, name)), exist_ok=True)
        with open(os.path.join(staging, name), "w", encoding="utf-8") as f:
            f.write(content)
    for name, arr in tiffs.items():
        os.makedirs(os.path.dirname(os.path.join(staging, name)), exist_ok=True)
        ds = gdal.GetDriverByName("GTiff").Create(os.path.join(staging, name), SAMPLES, LINES, 1, gdal.GDT_UInt16)
        ds.GetRasterBand(1).WriteArray(arr)
        ds = None
    if not as_zip:
        return os.path.join(staging, SAFE), dn
    zip_path = os.path.join(folder, SAFE + ".zip")
    with zipfile.ZipFile(zip_path, "w") as z:
        for root, _, names in os.walk(staging):
            for n in names:
                full = os.path.join(root, n)
                z.write(full, os.path.relpath(full, staging).replace(os.sep, "/"))
    return zip_path, dn


@unittest.skipUnless(HAVE_GDAL, "needs GDAL")
class TestSafeProduct(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from sar_flood_detection.safe import SafeProduct
        cls.SafeProduct = SafeProduct
        cls.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.zip_path, cls.dn = build_safe(cls.tmp.name)
        cls.folder_path, _ = build_safe(os.path.join(cls.tmp.name, "unzipped"), as_zip=False)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_finds_both_polarisations_in_zip_and_folder(self):
        for path in (self.zip_path, self.folder_path, os.path.join(self.folder_path, "manifest.safe")):
            with self.subTest(path=path):
                self.assertEqual(self.SafeProduct(path).polarisations, ["VH", "VV"])

    def test_reads_window_from_zip(self):
        p = self.SafeProduct(self.zip_path)
        np.testing.assert_array_equal(p.read_dn("vh", (10, 5, 30, 20)), self.dn["VH"][5:25, 10:40])

    def test_info_and_window(self):
        p = self.SafeProduct(self.zip_path)
        self.assertEqual(p.info("VV")["samples"], SAMPLES)
        self.assertEqual(p.window_for_bbox("VV", (105.01, 19.97, 105.02, 19.98), margin=0), (10, 20, 10, 10))

    def test_backscatter_matches_formula(self):
        p = self.SafeProduct(self.zip_path)
        out = p.backscatter("VH", (0, 0, SAMPLES, 1))[0]
        cols = np.arange(SAMPLES)
        a = 500.0 + 100.0 * cols / SAMPLES
        noise = 100.0 - 50.0 * cols / SAMPLES
        expected = np.maximum((self.dn["VH"][0].astype(float) ** 2 - noise) / a ** 2, 1e-5)
        expected[:3] = np.nan
        np.testing.assert_allclose(out, expected, rtol=1e-5)

    def test_strips_give_the_same_result(self):
        p = self.SafeProduct(self.zip_path)
        window = (2, 3, SAMPLES - 4, 17)
        np.testing.assert_array_equal(p.backscatter("VH", window, strip_rows=4),
                                      p.backscatter("VH", window, strip_rows=1000))

    def test_missing_polarisation_raises(self):
        with self.assertRaises(ValueError):
            self.SafeProduct(self.zip_path).read_dn("HH")

    def test_not_a_product_raises(self):
        empty = os.path.join(self.tmp.name, "empty")
        os.makedirs(empty, exist_ok=True)
        with self.assertRaises(IOError):
            self.SafeProduct(empty)


if __name__ == "__main__":
    unittest.main()
