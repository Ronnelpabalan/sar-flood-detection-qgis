"""End-to-end tests of the Processing tools on synthetic GeoTIFFs.

These need QGIS's Python (e.g. python-qgis-ltr.bat) and are skipped elsewhere.
"""

import os
import sys
import tempfile
import unittest

import numpy as np

try:
    from osgeo import gdal, osr
    from qgis.core import QgsApplication, QgsProcessingException, QgsProcessingFeedback
    HAVE_QGIS = True
except ImportError:
    HAVE_QGIS = False

ENL = 4.4
GEOTRANSFORM = (-1.15, 9e-5, 0.0, 9.47, 0.0, -9e-5)     # ~10 m pixels in EPSG:4326, like Sen1Floods11


def speckled(true_db, rng):
    linear = 10 ** (true_db / 10) * rng.gamma(ENL, 1 / ENL, true_db.shape)
    return (10 * np.log10(linear)).astype("float32")


def write_raster(path, bands, geotransform=GEOTRANSFORM, nodata=None):
    height, width = bands[0].shape
    ds = gdal.GetDriverByName("GTiff").Create(path, width, height, len(bands), gdal.GDT_Float32)
    ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    for i, arr in enumerate(bands, start=1):
        b = ds.GetRasterBand(i)
        b.WriteArray(arr)
        if nodata is not None:
            b.SetNoDataValue(nodata)
    ds = None
    return path


def read_output(path):
    ds = gdal.Open(path)
    return ds.GetRasterBand(1).ReadAsArray(), ds.GetGeoTransform(), ds.GetRasterBand(1).GetNoDataValue()


if HAVE_QGIS:
    class RecordingFeedback(QgsProcessingFeedback):
        """Feedback that remembers warnings, so tests can check them."""

        def __init__(self):
            super().__init__()
            self.warnings = []

        def pushWarning(self, warning):
            self.warnings.append(warning)
            super().pushWarning(warning)


@unittest.skipUnless(HAVE_QGIS, "needs QGIS Python")
class TestAlgorithms(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QgsApplication([], False)
        cls.app.initQgis()
        sys.path.append(os.path.join(QgsApplication.prefixPath(), "python", "plugins"))
        from processing.core.Processing import Processing
        import processing
        Processing.initialize()
        cls.run_alg = staticmethod(processing.run)

        import sar_flood_detection
        cls.plugin = sar_flood_detection.classFactory(None)
        cls.plugin.initProcessing()

        cls.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)   # QGIS may keep inputs open on Windows
        rng = np.random.default_rng(0)
        # 120 x 120 chip: land everywhere, a water square, and a block with no data
        vv_true = np.full((120, 120), -10.0)
        vh_true = np.full((120, 120), -17.0)
        vv_true[30:70, 30:70] = -21.0
        vh_true[30:70, 30:70] = -28.0
        cls.water = np.zeros((120, 120), bool)
        cls.water[30:70, 30:70] = True
        vv, vh = speckled(vv_true, rng), speckled(vh_true, rng)
        vv[100:110, 100:110] = np.nan
        vh[100:110, 100:110] = np.nan
        cls.db_path = write_raster(os.path.join(cls.tmp.name, "chip_db.tif"), [vv, vh])
        cls.linear_path = write_raster(os.path.join(cls.tmp.name, "chip_linear.tif"),
                                       [10 ** (vv / 10), 10 ** (vh / 10)])
        cls.land_only_path = write_raster(os.path.join(cls.tmp.name, "land_only.tif"),
                                          [speckled(np.full((120, 120), -17.0), rng)])

    @classmethod
    def tearDownClass(cls):
        cls.plugin.unload()
        cls.tmp.cleanup()

    def out(self, name):
        return os.path.join(self.tmp.name, name)

    # ---- water map ----

    def water_map(self, output, **params):
        defaults = {"INPUT": self.db_path, "BAND": 2, "UNITS": 0, "FILTER": 1, "WINDOW": 5,
                    "METHOD": 0, "THRESHOLD": -23.3, "OUTPUT": self.out(output)}
        defaults.update(params)
        return self.run_alg("sarflood:watermapsingleimage", defaults)

    def test_fixed_threshold_finds_the_water_square(self):
        result = self.water_map("fixed.tif")
        mask, gt, nodata = read_output(result["OUTPUT"])
        self.assertEqual(nodata, 255)
        self.assertEqual(gt, GEOTRANSFORM)
        self.assertTrue((mask[100:110, 100:110] == 255).all())
        self.assertGreater(mask[self.water].mean(), 0.97)
        land = ~self.water & (mask != 255)
        self.assertLess(mask[land].mean(), 0.02)
        self.assertEqual(result["THRESHOLD_USED"], -23.3)

    def test_linear_input_gives_same_mask(self):
        a = read_output(self.water_map("from_db.tif")["OUTPUT"])[0]
        b = read_output(self.water_map("from_linear.tif", INPUT=self.linear_path, UNITS=1)["OUTPUT"])[0]
        self.assertGreater(np.mean(a == b), 0.999)

    def test_each_filter_option_runs(self):
        for f in (0, 1, 2):
            with self.subTest(filter=f):
                mask = read_output(self.water_map(f"filter_{f}.tif", FILTER=f)["OUTPUT"])[0]
                self.assertGreater(mask[self.water].mean(), 0.9)

    def test_otsu_on_water_and_land(self):
        fb = RecordingFeedback()
        params = {"INPUT": self.db_path, "BAND": 2, "UNITS": 0, "FILTER": 1, "WINDOW": 5,
                  "METHOD": 1, "THRESHOLD": 0, "OUTPUT": self.out("otsu.tif")}
        result = self.run_alg("sarflood:watermapsingleimage", params, feedback=fb)
        self.assertTrue(-28 < result["THRESHOLD_USED"] < -17)
        self.assertEqual(fb.warnings, [])

    def test_otsu_warns_without_water(self):
        fb = RecordingFeedback()
        params = {"INPUT": self.land_only_path, "BAND": 1, "UNITS": 0, "FILTER": 1, "WINDOW": 5,
                  "METHOD": 1, "THRESHOLD": 0, "OUTPUT": self.out("otsu_land.tif")}
        self.run_alg("sarflood:watermapsingleimage", params, feedback=fb)
        self.assertEqual(len(fb.warnings), 1)

    def test_even_window_is_rejected(self):
        with self.assertRaises(QgsProcessingException):
            self.water_map("even.tif", WINDOW=4)

    # ---- change detection ----

    def test_change_detection_classes(self):
        rng = np.random.default_rng(1)
        pre_true = np.full((120, 120), -17.0)
        pre_true[:, 0:15] = -28.0                    # river, present before and after
        post_true = pre_true.copy()
        post_true[40:80, 40:80] = -28.0              # newly flooded square
        pre = write_raster(self.out("pre.tif"), [speckled(pre_true, rng)])
        post = write_raster(self.out("post.tif"), [speckled(post_true, rng)])
        result = self.run_alg("sarflood:floodchangedetection", {
            "PRE": pre, "PRE_BAND": 1, "POST": post, "POST_BAND": 1, "UNITS": 0, "FILTER": 1,
            "WINDOW": 5, "WATER_THRESHOLD": -23.3, "MIN_DROP": 3.0, "OUTPUT": self.out("change.tif")})
        classes = read_output(result["OUTPUT"])[0]
        self.assertGreater(np.mean(classes[42:78, 42:78] == 1), 0.97)    # flood (inner part)
        self.assertGreater(np.mean(classes[:, 2:13] == 2), 0.97)         # permanent water
        self.assertGreater(np.mean(classes[90:, 30:] == 0), 0.97)        # dry
        self.assertAlmostEqual(result["FLOOD_PERCENT"], 100 * 40 * 40 / 120 ** 2, delta=1.0)

    def test_change_detection_rejects_different_grids(self):
        rng = np.random.default_rng(2)
        pre = write_raster(self.out("pre_a.tif"), [speckled(np.full((50, 50), -17.0), rng)])
        shifted = (GEOTRANSFORM[0] + 0.01,) + GEOTRANSFORM[1:]
        post = write_raster(self.out("post_b.tif"), [speckled(np.full((50, 50), -17.0), rng)], shifted)
        with self.assertRaises(QgsProcessingException):
            self.run_alg("sarflood:floodchangedetection", {
                "PRE": pre, "PRE_BAND": 1, "POST": post, "POST_BAND": 1,
                "OUTPUT": self.out("bad.tif")})


if __name__ == "__main__":
    unittest.main()
