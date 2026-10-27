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

    # ---- flood confidence ----

    def confidence_inputs(self):
        """Pre/post pair with a river (columns 0-14), a flooded square on flat land (rows and columns
        40-79) and a darkened patch on a steep hill (rows 90-109, columns 90-109). Written once: QGIS
        keeps input files open on Windows, so they cannot be rewritten."""
        if hasattr(type(self), "_confidence_inputs"):
            return type(self)._confidence_inputs
        rng = np.random.default_rng(3)
        pre_true = np.full((120, 120), -17.0)
        pre_true[:, 0:15] = -28.0
        post_true = pre_true.copy()
        post_true[40:80, 40:80] = -28.0
        post_true[90:110, 90:110] = -28.0
        height = np.zeros((120, 120), dtype="float32")
        height[:, 85:] = (np.arange(35, dtype="float32") * 10.0)[None, :]    # 45 degree slope, up to 340 m
        mask = np.zeros((120, 120), dtype="float32")
        mask[45:50, 45:50] = 2                                                # shadow inside the flood
        type(self)._confidence_inputs = {
            "PRE": write_raster(self.out("c_pre.tif"), [speckled(pre_true, rng)]),
            "POST": write_raster(self.out("c_post.tif"), [speckled(post_true, rng)]),
            "dem": write_raster(self.out("c_dem.tif"), [height]),
            "mask": write_raster(self.out("c_mask.tif"), [mask]),
        }
        return type(self)._confidence_inputs

    def confidence(self, inputs, name, **params):
        defaults = {"PRE": inputs["PRE"], "PRE_BAND": 1, "POST": inputs["POST"], "POST_BAND": 1,
                    "UNITS": 0, "FILTER": 1, "WINDOW": 5, "WATER_THRESHOLD": -23.3, "MIN_DROP": 3.0,
                    "OUTPUT": self.out(f"{name}.tif"), "SCORE": self.out(f"{name}_score.tif")}
        defaults.update(params)
        return self.run_alg("sarflood:floodconfidence", defaults)

    def test_confidence_without_dem_matches_change_detection(self):
        inputs = self.confidence_inputs()
        result = self.confidence(inputs, "conf_plain")
        change = self.run_alg("sarflood:floodchangedetection", {
            "PRE": inputs["PRE"], "PRE_BAND": 1, "POST": inputs["POST"], "POST_BAND": 1, "UNITS": 0,
            "FILTER": 1, "WINDOW": 5, "WATER_THRESHOLD": -23.3, "MIN_DROP": 3.0,
            "OUTPUT": self.out("conf_change.tif")})
        classes, old = read_output(result["OUTPUT"])[0], read_output(change["OUTPUT"])[0]
        np.testing.assert_array_equal(np.isin(classes, (2, 3)), old == 1)
        np.testing.assert_array_equal(classes == 4, old == 2)
        score = read_output(result["SCORE"])[0]
        self.assertTrue(np.nanmin(score) >= 0 and np.nanmax(score) <= 1)
        self.assertGreater(np.mean(classes[42:78, 42:78] == 3), 0.97)     # flood square: high

    def test_confidence_with_dem_and_mask(self):
        inputs = self.confidence_inputs()
        result = self.confidence(inputs, "conf_terrain", DEM=[inputs["dem"]], POST_MASK=inputs["mask"])
        classes = read_output(result["OUTPUT"])[0]
        self.assertTrue((classes[45:50, 45:50] == 255).all())             # shadow: no data
        self.assertGreater(np.mean(classes[52:78, 52:78] == 3), 0.97)     # flat flood: still high
        self.assertGreater(np.mean(classes[92:108, 92:108] == 0), 0.97)   # hill: not flood
        plain = read_output(self.confidence(inputs, "conf_plain2")["OUTPUT"])[0]
        self.assertGreater(np.mean(plain[92:108, 92:108] == 3), 0.97)     # ...which it was without the DEM
        self.assertAlmostEqual(result["PERMANENT_KM2"], 120 * 15 * 0.01 * 0.0099, delta=0.02)
        alg = QgsApplication.processingRegistry().algorithmById("sarflood:floodconfidence")
        self.assertEqual(alg.parameterDefinition("SLOPE_WINDOW").defaultValue(), 90.0)
        self.assertEqual(alg.parameterDefinition("NEAR_WATER").defaultValue(), 3.0)

    def test_confidence_rejects_mask_on_other_grid(self):
        inputs = self.confidence_inputs()
        shifted = (GEOTRANSFORM[0] + 0.01,) + GEOTRANSFORM[1:]
        bad = write_raster(self.out("c_mask_shifted.tif"), [np.zeros((120, 120), "float32")], shifted)
        with self.assertRaises(QgsProcessingException):
            self.confidence(inputs, "conf_bad", PRE_MASK=bad)

    # ---- prepare GRD (the full run on real data is checked by hand; see README) ----

    def test_prepare_grd_is_registered(self):
        alg = QgsApplication.processingRegistry().algorithmById("sarflood:preparegrd")
        self.assertIsNotNone(alg)
        names = {p.name() for p in alg.parameterDefinitions()}
        self.assertTrue({"INPUT", "EXTENT", "POLARISATIONS", "SPACING", "ORBIT", "DATA_FOLDER",
                         "DOWNLOAD", "OUTPUT", "MASK"} <= names)

    def test_prepare_grd_rejects_other_files(self):
        with self.assertRaises(QgsProcessingException):
            self.run_alg("sarflood:preparegrd", {
                "INPUT": self.db_path, "EXTENT": "105.4,106.0,19.6,20.2 [EPSG:4326]",
                "DOWNLOAD": False, "OUTPUT": self.out("not_grd.tif")})

    # ---- styling when QGIS loads the result ----

    def test_loaded_output_has_transparent_not_water_class(self):
        from qgis.core import (QgsProcessingContext, QgsProcessingOutputLayerDefinition,
                               QgsProject, QgsRasterLayer)
        context = QgsProcessingContext()
        context.setProject(QgsProject.instance())
        path = self.out("styled.tif")
        params = {"INPUT": self.db_path, "BAND": 2, "UNITS": 0, "FILTER": 1, "WINDOW": 5,
                  "METHOD": 0, "THRESHOLD": -23.3,
                  "OUTPUT": QgsProcessingOutputLayerDefinition(path, QgsProject.instance())}
        alg = QgsApplication.processingRegistry().createAlgorithmById("sarflood:watermapsingleimage")
        results, ok = alg.run(params, context, QgsProcessingFeedback())
        self.assertTrue(ok)
        self.assertTrue(context.willLoadLayerOnCompletion(results["OUTPUT"]))
        styler = context.layerToLoadOnCompletionDetails(results["OUTPUT"]).postProcessor()
        self.assertIsNotNone(styler)

        layer = QgsRasterLayer(results["OUTPUT"], "Water mask")
        styler.postProcessLayer(layer, context, QgsProcessingFeedback())
        renderer = layer.renderer()
        self.assertEqual(renderer.type(), "paletted")
        classes = {c.value: c for c in renderer.classes()}
        self.assertEqual(classes[0].label, "Not water")
        self.assertEqual(classes[0].color.alpha(), 0)          # transparent over the SAR image
        self.assertEqual(classes[1].label, "Water")
        self.assertEqual(classes[1].color.alpha(), 255)

    def test_gray_styler_shows_one_band_stretched(self):
        from qgis.core import QgsProcessingContext, QgsRasterLayer
        from sar_flood_detection.algorithms.common import GrayStyler
        layer = QgsRasterLayer(self.db_path, "Backscatter")
        GrayStyler(2).postProcessLayer(layer, QgsProcessingContext(), QgsProcessingFeedback())
        renderer = layer.renderer()
        self.assertEqual(renderer.type(), "singlebandgray")
        self.assertEqual(renderer.inputBand(), 2)
        ce = renderer.contrastEnhancement()
        ds = gdal.Open(self.db_path)              # keep a reference, or GDAL closes the file early
        vh = ds.GetRasterBand(2).ReadAsArray()
        low, high = np.nanpercentile(vh, [2, 98])
        self.assertAlmostEqual(ce.minimumValue(), low, delta=1.0)     # QGIS estimates from a histogram
        self.assertAlmostEqual(ce.maximumValue(), high, delta=1.0)

    def test_score_styler_uses_a_colour_ramp(self):
        from qgis.core import QgsProcessingContext, QgsRasterLayer
        from sar_flood_detection.algorithms.common import ScoreStyler
        from sar_flood_detection.algorithms.flood_confidence import SCORE_COLORS
        layer = QgsRasterLayer(self.db_path, "Score")
        ScoreStyler(SCORE_COLORS).postProcessLayer(layer, QgsProcessingContext(), QgsProcessingFeedback())
        renderer = layer.renderer()
        self.assertEqual(renderer.type(), "singlebandpseudocolor")
        self.assertEqual((renderer.classificationMin(), renderer.classificationMax()), (0.0, 1.0))
        items = renderer.shader().rasterShaderFunction().colorRampItemList()
        self.assertEqual(len(items), len(SCORE_COLORS))
        self.assertEqual(items[0].color.alpha(), 0)            # score 0 is transparent

    def test_output_file_palette_is_not_black(self):
        mask_path = self.water_map("palette.tif")["OUTPUT"]
        ds = gdal.Open(mask_path)                 # keep a reference, or GDAL closes the file early
        table = ds.GetRasterBand(1).GetRasterColorTable()
        self.assertEqual(table.GetColorEntry(0)[:3], (255, 255, 255))   # not water: white in the file


if __name__ == "__main__":
    unittest.main()
