"""Sentinel-1 orbit state vectors: reading and interpolation.

Two sources give the satellite position and velocity (Earth-fixed, metres):

- the product annotation (generalAnnotation/orbitList): 10 s steps around
  the image, produced when the image was processed;
- precise orbit files (AUX_POEORB, .EOF), published by ESA about 20 days
  later; restituted files (AUX_RESORB) are a quicker, less accurate version.

Between the 10 s steps the position is interpolated with a cubic Hermite
spline, which uses both the positions and the velocities.

Times are seconds since 1970-01-01 UTC (float).
"""

import re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

import numpy as np
from scipy.interpolate import CubicHermiteSpline

# S1D_OPER_AUX_POEORB_OPOD_20261023T071430_V20261002T225942_20261004T005942.EOF(.zip)
_EOF_NAME = re.compile(
    r"(?P<mission>S1[A-D])_OPER_AUX_(?P<kind>POEORB|RESORB)_OPOD_(?P<produced>\d{8}T\d{6})_"
    r"V(?P<start>\d{8}T\d{6})_(?P<stop>\d{8}T\d{6})\.EOF(\.zip)?$")


def parse_time(text):
    """UTC time text ('2026-10-03T22:50:43.258500', optionally 'UTC=...') to seconds since 1970."""
    text = text.strip()
    if text.startswith("UTC="):
        text = text[4:]
    fmt = "%Y%m%dT%H%M%S" if "-" not in text else ("%Y-%m-%dT%H:%M:%S.%f" if "." in text else "%Y-%m-%dT%H:%M:%S")
    return datetime.strptime(text, fmt).replace(tzinfo=timezone.utc).timestamp()


def _state_vectors(times, pos, vel):
    order = np.argsort(times)
    return {"times": np.array(times)[order], "positions": np.array(pos)[order],
            "velocities": np.array(vel)[order]}


def parse_eof(xml_text):
    """State vectors from a precise or restituted orbit file (.EOF)."""
    root = ET.fromstring(xml_text)
    times, pos, vel = [], [], []
    for osv in root.iter("OSV"):
        times.append(parse_time(osv.findtext("UTC")))
        pos.append([float(osv.findtext(k)) for k in ("X", "Y", "Z")])
        vel.append([float(osv.findtext(k)) for k in ("VX", "VY", "VZ")])
    if not times:
        raise ValueError("No OSV elements found in the orbit file.")
    return _state_vectors(times, pos, vel)


def parse_annotation_orbit(xml_text):
    """State vectors from a product annotation file (annotation/s1*.xml)."""
    root = ET.fromstring(xml_text)
    times, pos, vel = [], [], []
    for o in root.iter("orbit"):
        times.append(parse_time(o.findtext("time")))
        pos.append([float(o.findtext(f"position/{k}")) for k in "xyz"])
        vel.append([float(o.findtext(f"velocity/{k}")) for k in "xyz"])
    if not times:
        raise ValueError("No orbit elements found in the annotation.")
    return _state_vectors(times, pos, vel)


def parse_eof_name(name):
    """Mission, kind, production time and validity period from an orbit file name, or None."""
    m = _EOF_NAME.search(name)
    if not m:
        return None
    return {"name": name, "mission": m.group("mission"), "kind": m.group("kind"),
            "produced": parse_time(m.group("produced")),
            "start": parse_time(m.group("start")), "stop": parse_time(m.group("stop"))}


def select_orbit_file(names, mission, start, stop, margin=60.0):
    """Best orbit file name that covers [start - margin, stop + margin], or None.

    Precise (POEORB) files are preferred over restituted (RESORB); among
    equals, the most recently produced one.
    """
    best = None
    for n in names:
        f = parse_eof_name(n)
        if f is None or f["mission"] != mission:
            continue
        if f["start"] > start - margin or f["stop"] < stop + margin:
            continue
        key = (f["kind"] == "POEORB", f["produced"])
        if best is None or key > best[0]:
            best = (key, n)
    return None if best is None else best[1]


class Orbit:
    """Interpolated satellite position and velocity.

    Only state vectors within margin seconds of [start, stop] are kept, so
    a 26-hour precise orbit file is reduced to the part around the image.
    """

    def __init__(self, vectors, start=None, stop=None, margin=60.0):
        t = vectors["times"]
        keep = np.ones(len(t), bool)
        if start is not None:
            keep &= t >= start - margin
        if stop is not None:
            keep &= t <= stop + margin
        if keep.sum() < 4:
            raise ValueError("The orbit data do not cover the image time.")
        self.times = t[keep]
        self._x = CubicHermiteSpline(self.times, vectors["positions"][keep], vectors["velocities"][keep], axis=0)
        self._v = self._x.derivative()
        self._a = self._x.derivative(2)

    def _check(self, t):
        t = np.asarray(t, dtype="float64")
        if t.min() < self.times[0] or t.max() > self.times[-1]:
            raise ValueError("Time outside the orbit data.")
        return t

    def position(self, t):
        """Earth-fixed position (m) at time(s) t; shape (..., 3)."""
        return self._x(self._check(t))

    def velocity(self, t):
        """Earth-fixed velocity (m/s) at time(s) t; shape (..., 3)."""
        return self._v(self._check(t))

    def acceleration(self, t):
        """Earth-fixed acceleration (m/s²) at time(s) t; shape (..., 3)."""
        return self._a(self._check(t))
