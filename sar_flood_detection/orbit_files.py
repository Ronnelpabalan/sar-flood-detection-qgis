"""Finding, downloading and caching Sentinel-1 orbit files.

Files come from ESA's STEP server, the source SNAP also uses:
https://step.esa.int/auxdata/orbits/Sentinel-1/<POEORB|RESORB>/<S1x>/<year>/<month>/
Each month folder holds the files whose validity starts in that month.
"""

import os
import re
import urllib.request
import zipfile
from datetime import datetime, timedelta, timezone

from .core.orbit import parse_eof, parse_eof_name, select_orbit_file

BASE_URL = "https://step.esa.int/auxdata/orbits/Sentinel-1"
KINDS = ("POEORB", "RESORB")


def default_cache_dir():
    return os.path.join(os.path.expanduser("~"), ".sar_flood_detection", "orbits")


def _fetch(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def list_remote(kind, mission, year, month, fetch=_fetch):
    """Orbit file names in one month folder of the STEP server."""
    url = f"{BASE_URL}/{kind}/{mission}/{year}/{month:02d}/"
    html = fetch(url).decode("utf-8", "replace")
    return sorted(set(re.findall(r'href="(S1[A-D]_OPER_AUX_[A-Z]+_OPOD_[^"/]+\.EOF(?:\.zip)?)"', html)))


def find_orbit_file(mission, start, stop, cache_dir=None, download=True, fetch=_fetch, feedback=None):
    """Local path of the best orbit file for an image (start/stop: seconds since 1970 UTC).

    Looks in the cache first; otherwise lists the server folders and
    downloads the file into the cache. Returns None if no file is found.
    """
    cache_dir = cache_dir or default_cache_dir()
    if os.path.isdir(cache_dir):
        name = select_orbit_file(os.listdir(cache_dir), mission, start, stop)
        if name and name.split("_")[3] == "POEORB":
            return os.path.join(cache_dir, name)
    else:
        name = None
    if not download:
        return os.path.join(cache_dir, name) if name else None

    day = datetime.fromtimestamp(start, timezone.utc)
    months = sorted({(d.year, d.month) for d in (day - timedelta(days=1), day)})
    for kind in KINDS:
        names = []
        for year, month in months:
            try:
                names += list_remote(kind, mission, year, month, fetch)
            except OSError as e:
                if feedback is not None:
                    feedback.pushInfo(f"Cannot list {kind} {year}-{month:02d}: {e}")
        best = select_orbit_file(names, mission, start, stop)
        if best:
            valid_from = datetime.fromtimestamp(parse_eof_name(best)["start"], timezone.utc)
            url = f"{BASE_URL}/{kind}/{mission}/{valid_from.year}/{valid_from.month:02d}/{best}"
            if feedback is not None:
                feedback.pushInfo(f"Downloading orbit file {best}")
            data = fetch(url)
            os.makedirs(cache_dir, exist_ok=True)
            path = os.path.join(cache_dir, best)
            with open(path, "wb") as f:
                f.write(data)
            return path
    return os.path.join(cache_dir, name) if name else None


def read_orbit_file(path):
    """State vectors from a .EOF file or a .EOF.zip containing one."""
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            eof = [n for n in z.namelist() if n.upper().endswith(".EOF")]
            if not eof:
                raise IOError(f"No .EOF file inside {path}")
            text = z.read(eof[0]).decode("utf-8")
    else:
        with open(path, encoding="utf-8") as f:
            text = f.read()
    return parse_eof(text)
