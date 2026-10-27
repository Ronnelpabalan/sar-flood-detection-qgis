"""Tests for core/orbit.py and orbit_files.py (no QGIS or internet needed)."""

import os
import tempfile
import unittest
import zipfile
from datetime import datetime, timedelta, timezone

import numpy as np

from sar_flood_detection.core.orbit import (Orbit, parse_annotation_orbit, parse_eof, parse_eof_name,
                                            parse_time, select_orbit_file)
from sar_flood_detection.orbit_files import find_orbit_file, read_orbit_file
from tests.s1_fixtures import annotation_orbit_xml, eof_xml

T0 = datetime(2026, 10, 3, 22, 49, 0, tzinfo=timezone.utc)
RADIUS, PERIOD, INCLINATION = 7.07e6, 5924.0, np.radians(98.2)


def circular(t):
    """Exact position and velocity on an inclined circular orbit (t: seconds after T0)."""
    w = 2 * np.pi / PERIOD
    a = w * np.asarray(t, dtype="float64")
    p = RADIUS * np.stack([np.cos(a), np.sin(a) * np.cos(INCLINATION), np.sin(a) * np.sin(INCLINATION)], -1)
    v = RADIUS * w * np.stack([-np.sin(a), np.cos(a) * np.cos(INCLINATION), np.cos(a) * np.sin(INCLINATION)], -1)
    return p, v


def utc(seconds_after_t0):
    return (T0 + timedelta(seconds=float(seconds_after_t0))).strftime("%Y-%m-%dT%H:%M:%S.%f")


def name(kind, produced, start, stop, mission="S1D"):
    return f"{mission}_OPER_AUX_{kind}_OPOD_{produced}_V{start}_{stop}.EOF.zip"


class TestOrbit(unittest.TestCase):
    def setUp(self):
        steps = np.arange(0, 301, 10.0)                       # 10 s state vectors, like ESA's files
        p, v = circular(steps)
        self.vectors = parse_eof(eof_xml([utc(s) for s in steps], p, v))
        self.t0 = T0.timestamp()

    def test_parse_time(self):
        self.assertEqual(parse_time("UTC=2026-10-03T22:49:00.000000"), self.t0)
        self.assertEqual(parse_time("2026-10-03T22:49:00"), self.t0)
        self.assertEqual(parse_time("20261003T224900"), self.t0)
        self.assertAlmostEqual(parse_time("2026-10-03T22:49:00.258500") - self.t0, 0.2585, places=6)

    def test_interpolation_matches_exact_orbit(self):
        orbit = Orbit(self.vectors)
        t = np.random.default_rng(0).uniform(5, 295, 200)
        p, v = circular(t)
        self.assertLess(np.abs(orbit.position(self.t0 + t) - p).max(), 1e-3)      # metres
        self.assertLess(np.abs(orbit.velocity(self.t0 + t) - v).max(), 1e-3)      # m/s

    def test_trims_to_image_time(self):
        orbit = Orbit(self.vectors, self.t0 + 100, self.t0 + 130, margin=40)
        self.assertEqual((orbit.times[0] - self.t0, orbit.times[-1] - self.t0), (60.0, 170.0))

    def test_time_outside_orbit_raises(self):
        with self.assertRaises(ValueError):
            Orbit(self.vectors).position(self.t0 + 400)
        with self.assertRaises(ValueError):
            Orbit(self.vectors, self.t0 + 1000, self.t0 + 1030)

    def test_annotation_orbit_matches_eof(self):
        steps = np.arange(0, 151, 10.0)
        p, v = circular(steps)
        ann = parse_annotation_orbit(annotation_orbit_xml([utc(s) for s in steps], p, v))
        np.testing.assert_allclose(ann["positions"], self.vectors["positions"][:16], atol=1e-6)


class TestOrbitFileChoice(unittest.TestCase):
    start = parse_time("2026-10-03T22:50:43")
    stop = parse_time("2026-10-03T22:51:08")

    def test_parse_name(self):
        f = parse_eof_name(name("POEORB", "20261023T071430", "20261002T225942", "20261004T005942"))
        self.assertEqual((f["mission"], f["kind"]), ("S1D", "POEORB"))
        self.assertEqual(f["start"], parse_time("2026-10-02T22:59:42"))
        self.assertIsNone(parse_eof_name("readme.txt"))

    def test_picks_file_covering_image(self):
        names = [name("POEORB", "20261024T071421", "20261003T225942", "20261005T005942"),   # starts after image
                 name("POEORB", "20261023T071430", "20261002T225942", "20261004T005942"),
                 name("POEORB", "20261023T071430", "20261002T225942", "20261004T005942", "S1A")]
        self.assertEqual(select_orbit_file(names, "S1D", self.start, self.stop), names[1])

    def test_prefers_precise_then_latest(self):
        res = name("RESORB", "20261004T010000", "20261003T220000", "20261004T010000")
        old = name("POEORB", "20261022T071440", "20261002T225942", "20261004T005942")
        new = name("POEORB", "20261023T071430", "20261002T225942", "20261004T005942")
        self.assertEqual(select_orbit_file([res, old, new], "S1D", self.start, self.stop), new)
        self.assertEqual(select_orbit_file([res], "S1D", self.start, self.stop), res)
        self.assertIsNone(select_orbit_file([], "S1D", self.start, self.stop))


class FakeServer:
    """Stands in for step.esa.int: month folders hold files by start of validity."""

    def __init__(self, files):
        self.files, self.requests = files, []

    def __call__(self, url):
        self.requests.append(url)
        folder, _, filename = url.rpartition("/")
        if filename:
            return self.files[url]
        names = [u.rpartition("/")[2] for u in self.files if u.rpartition("/")[0] == folder]
        return "".join(f'<a href="{n}">{n}</a>' for n in names).encode()


class TestFindOrbitFile(unittest.TestCase):
    start = parse_time("2026-10-01T22:50:43")
    stop = parse_time("2026-10-01T22:51:08")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        base = "https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/S1D/2026"
        self.wanted = name("POEORB", "20261021T071430", "20260930T225942", "20261002T005942")
        self.server = FakeServer({
            f"{base}/09/{self.wanted}": b"zip-bytes",
            f"{base}/10/" + name("POEORB", "20261022T071440", "20261001T225942", "20261003T005942"): b"x",
        })

    def tearDown(self):
        self.tmp.cleanup()

    def test_downloads_from_validity_month_and_caches(self):
        path = find_orbit_file("S1D", self.start, self.stop, self.tmp.name, fetch=self.server)
        self.assertEqual(os.path.basename(path), self.wanted)
        self.assertTrue(self.server.requests[-1].endswith(f"/2026/09/{self.wanted}"))
        with open(path, "rb") as f:
            self.assertEqual(f.read(), b"zip-bytes")

        def no_network(url):
            raise AssertionError("should use the cache")
        self.assertEqual(find_orbit_file("S1D", self.start, self.stop, self.tmp.name, fetch=no_network), path)

    def test_nothing_found_returns_none(self):
        late = parse_time("2026-12-01T22:50:43")
        def empty(url):
            return b"<html></html>"
        self.assertIsNone(find_orbit_file("S1D", late, late + 25, self.tmp.name, fetch=empty))

    def test_read_zipped_eof(self):
        steps = np.arange(0, 41, 10.0)
        p, v = circular(steps)
        path = os.path.join(self.tmp.name, "orbit.EOF.zip")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("S1D_OPER_AUX_POEORB.EOF", eof_xml([utc(s) for s in steps], p, v))
        self.assertEqual(read_orbit_file(path)["positions"].shape, (5, 3))


if __name__ == "__main__":
    unittest.main()
