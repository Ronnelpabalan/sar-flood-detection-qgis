"""Tests for core/terrain.py and dem.py (GDAL from QGIS Python, no internet needed)."""

import os
import tempfile
import unittest
import urllib.error

import numpy as np
from osgeo import gdal, osr

from sar_flood_detection.core.terrain import LAYOVER, SHADOW, bilinear, layover_shadow
from sar_flood_detection.dem import (GridLonLat, find_dem_tiles, make_grid, tile_name, tile_names,
                                     warp_to_grid)

THANH_HOA = (105.4, 19.6, 106.0, 20.2)


class TestBilinear(unittest.TestCase):
    def setUp(self):
        r, c = np.mgrid[0:5, 0:6]
        self.image = (10.0 * r + c).astype(np.float32)              # a plane: bilinear is exact

    def test_exact_on_a_plane(self):
        rows, cols = np.array([0.0, 1.25, 3.5, 4.0]), np.array([0.0, 2.5, 4.75, 5.0])
        np.testing.assert_allclose(bilinear(self.image, rows, cols), 10 * rows + cols, atol=1e-5)

    def test_outside_and_nan(self):
        out = bilinear(self.image, np.array([-0.1, 2.0, 4.01, np.nan]), np.array([1.0, 6.0, 1.0, 1.0]))
        self.assertTrue(np.isnan(out).all())
        self.image[2, 2] = np.nan
        self.assertTrue(np.isnan(bilinear(self.image, np.array([1.5]), np.array([1.5])))[0])
        self.assertFalse(np.isnan(bilinear(self.image, np.array([0.5]), np.array([0.5])))[0])


class TestLayoverShadow(unittest.TestCase):
    def test_one_line(self):
        # points ordered by distance from the track
        distance = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        look = np.array([0.50, 0.60, 0.70, 0.65, 0.80])             # point 3 hides behind point 2
        rng = np.array([100.0, 110.0, 105.0, 120.0, 130.0])         # points 1 and 2 overlap
        flags = layover_shadow(np.zeros(5, int), distance, look, rng)
        np.testing.assert_array_equal(flags, [0, LAYOVER, LAYOVER, SHADOW, 0])

    def test_lines_are_independent_and_order_does_not_matter(self):
        group = np.array([1, 1, 0, 0])
        distance = np.array([2.0, 1.0, 1.0, 2.0])
        look = np.array([0.6, 0.5, 0.9, 0.7])                       # line 0: second point in shadow
        rng = np.array([110.0, 100.0, 90.0, 95.0])                  # would overlap only across lines
        np.testing.assert_array_equal(layover_shadow(group, distance, look, rng), [0, 0, 0, SHADOW])

    def test_separation_ignores_close_neighbours(self):
        # flat ground with a 3 m bump: points 2 and 3 are only 1 m apart, so the bump
        # flips their slant ranges, but they are closer than the 30 m DEM cell
        distance = np.array([0.0, 30.0, 60.0, 61.0, 90.0])
        look = np.array([0.50, 0.51, 0.52, 0.5201, 0.53])
        rng = np.array([100.0, 116.0, 132.0, 131.0, 148.0])
        self.assertTrue((layover_shadow(np.zeros(5, int), distance, look, rng) == LAYOVER).sum() == 2)
        np.testing.assert_array_equal(layover_shadow(np.zeros(5, int), distance, look, rng, separation=30.0), 0)
        # a real steep slope over more than 30 m is still found
        rng[4] = 120.0
        flags = layover_shadow(np.zeros(5, int), distance, look, rng, separation=30.0)
        np.testing.assert_array_equal(flags, [0, 0, LAYOVER, 0, LAYOVER])     # 132 m at 60 vs 120 m at 90

    def test_empty(self):
        self.assertEqual(layover_shadow(np.zeros(0, int), np.zeros(0), np.zeros(0), np.zeros(0)).size, 0)


class FakeLonLat:
    """A 'grid' whose longitude/latitude are linear in column/row."""
    grid = {"size": (200, 150)}

    def at(self, rows, cols):
        cc, rr = np.meshgrid(np.asarray(cols, float), np.asarray(rows, float))
        return 105.0 + cc * 1e-4, 20.0 - rr * 1e-4


class FakeGeometry:
    """Smooth stand-ins for the Range-Doppler results, as functions of x, y, z."""

    def observe(self, xyz):
        xyz = np.atleast_2d(xyz)
        return {"line": xyz[:, 0] / 10.0, "pixel": np.linalg.norm(xyz, axis=1) / 10.0}


class TestGeometryLattice(unittest.TestCase):
    def test_close_to_direct_evaluation(self):
        from sar_flood_detection.core.geometry import geodetic_to_ecef
        from sar_flood_detection.terrain_correction import GeometryLattice
        geo, lonlat = FakeGeometry(), FakeLonLat()
        lattice = GeometryLattice(geo, lonlat, 0.0, 600.0, step=16)
        rows, cols = np.arange(0, 150, 7), np.arange(0, 200, 3)
        heights = np.random.default_rng(2).uniform(0, 600, (len(rows), len(cols)))
        got = lattice(rows, cols, heights)
        lon, lat = lonlat.at(rows, cols)
        exact = geo.observe(geodetic_to_ecef(lat.ravel(), lon.ravel(), heights.ravel()))
        for key in ("line", "pixel"):
            self.assertEqual(got[key].shape, heights.shape)
            np.testing.assert_allclose(got[key].ravel(), exact[key], atol=0.01)   # 10 cm in metres / 10


class TestDemTiles(unittest.TestCase):
    def test_names(self):
        self.assertEqual(tile_name(19, 105), "Copernicus_DSM_COG_10_N19_00_E105_00_DEM")
        self.assertEqual(tile_name(-1, -72), "Copernicus_DSM_COG_10_S01_00_W072_00_DEM")
        self.assertEqual(tile_names(THANH_HOA), [tile_name(19, 105), tile_name(20, 105)])
        self.assertEqual(tile_names((105.2, 19.2, 105.8, 19.8)), [tile_name(19, 105)])

    def test_download_cache_and_sea_tiles(self):
        land, sea = tile_name(19, 105), tile_name(20, 105)
        requests = []

        def fetch(url):
            requests.append(url)
            if sea in url:
                raise urllib.error.HTTPError(url, 404, "Not Found", None, None)
            return b"tif-bytes"

        with tempfile.TemporaryDirectory() as tmp:
            paths = find_dem_tiles(THANH_HOA, tmp, fetch=fetch)
            self.assertEqual([os.path.basename(p) for p in paths], [land + ".tif"])
            self.assertEqual(requests[0], f"https://copernicus-dem-30m.s3.amazonaws.com/{land}/{land}.tif")

            def no_network(url):
                raise AssertionError("should use the cache")
            self.assertEqual(find_dem_tiles(THANH_HOA, tmp, fetch=no_network), paths)

    def test_missing_tiles_without_download_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(IOError):
                find_dem_tiles(THANH_HOA, tmp, download=False)


class TestGrid(unittest.TestCase):
    def test_snapped_and_repeatable(self):
        a, b = make_grid(THANH_HOA, 10.0), make_grid(THANH_HOA, 10.0)
        self.assertEqual(a, b)
        self.assertEqual(a["epsg"], 32648)                          # UTM 48N
        x0, dx, _, y1, _, dy = a["geotransform"]
        self.assertEqual((x0 % 10, y1 % 10, dx, dy), (0.0, 0.0, 10.0, -10.0))

    def test_lonlat_matches_exact_transform(self):
        grid = make_grid(THANH_HOA, 10.0)
        ll = GridLonLat(grid)
        lon, lat = ll.block(1000, 3)
        src = osr.SpatialReference()
        src.ImportFromEPSG(grid["epsg"])
        dst = osr.SpatialReference()
        dst.ImportFromEPSG(4326)
        for s in (src, dst):
            s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        t = osr.CoordinateTransformation(src, dst)
        x0, dx, _, y1, _, dy = grid["geotransform"]
        for row, col in ((1000, 0), (1001, 1234), (1002, grid["size"][0] - 1)):
            exact = t.TransformPoint(x0 + (col + 0.5) * dx, y1 + (row + 0.5) * dy)
            self.assertAlmostEqual(lon[row - 1000, col], exact[0], delta=1e-8)   # 1e-8 degree ~ 1 mm
            self.assertAlmostEqual(lat[row - 1000, col], exact[1], delta=1e-8)

    def test_warp_to_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dem.tif")
            ds = gdal.GetDriverByName("GTiff").Create(path, 20, 20, 1, gdal.GDT_Float32)
            ds.SetGeoTransform((105.0, 0.1, 0, 21.0, 0, -0.1))
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(4326)
            ds.SetProjection(srs.ExportToWkt())
            ds.GetRasterBand(1).Fill(5.0)
            ds = None
            grid = make_grid((105.5, 19.5, 105.6, 19.6), 100.0)
            np.testing.assert_allclose(warp_to_grid([path], grid), 5.0)
            self.assertTrue((warp_to_grid([], grid, fill=0.0) == 0).all())


if __name__ == "__main__":
    unittest.main()
