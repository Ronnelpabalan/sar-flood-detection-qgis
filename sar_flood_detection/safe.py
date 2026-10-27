"""Reading Sentinel-1 GRD products (.SAFE.zip or unzipped .SAFE folder) with GDAL."""

import os
import re
import zipfile

from osgeo import gdal

from .core.annotation import parse_geolocation_grid, parse_product_info, window_for_bbox
from .core.calibration import calibrate, parse_calibration, parse_noise

# e.g. measurement/s1a-iw-grd-vh-20220908t...-002.tiff, annotation/calibration/noise-s1a-iw-grd-vh-....xml
_FILE_PATTERN = re.compile(
    r"(?P<folder>measurement|annotation/calibration|annotation)/"
    r"(?P<prefix>calibration-|noise-)?s1[a-z]-[a-z0-9]+-grd-(?P<pol>[hv]{2})-.*\.(?P<ext>tiff|xml)$")


class SafeProduct:
    """A Sentinel-1 GRD product, read without unzipping.

    path: the .SAFE.zip file, the .SAFE folder, or its manifest.safe.
    """

    def __init__(self, path):
        path = os.path.abspath(path)
        if os.path.basename(path).lower() == "manifest.safe":
            path = os.path.dirname(path)
        self.path = path
        self.is_zip = zipfile.is_zipfile(path)
        if self.is_zip:
            with zipfile.ZipFile(path) as z:
                names = z.namelist()
        elif os.path.isdir(path):
            names = [os.path.relpath(os.path.join(root, f), path).replace(os.sep, "/")
                     for root, _, files in os.walk(path) for f in files]
        else:
            raise IOError(f"Not a Sentinel-1 .SAFE.zip or .SAFE folder: {path}")

        self.files = {}
        for name in names:
            m = _FILE_PATTERN.search(name)
            if not m:
                continue
            kind = {"measurement": "measurement", "annotation": "annotation"}.get(m.group("folder"))
            if kind is None:
                kind = m.group("prefix").rstrip("-")         # calibration or noise
            self.files.setdefault(m.group("pol").upper(), {})[kind] = name
        self.files = {p: f for p, f in self.files.items()
                      if {"measurement", "annotation", "calibration", "noise"} <= f.keys()}
        if not self.files:
            raise IOError(f"No GRD measurement with annotation, calibration and noise files in {path}")

    @property
    def polarisations(self):
        return sorted(self.files)

    def _check(self, pol):
        pol = pol.upper()
        if pol not in self.files:
            raise ValueError(f"Polarisation {pol} not in this product (has {', '.join(self.polarisations)}).")
        return pol

    def read_text(self, pol, kind):
        """Text of one XML file: kind is annotation, calibration or noise."""
        name = self.files[self._check(pol)][kind]
        if self.is_zip:
            with zipfile.ZipFile(self.path) as z:
                return z.read(name).decode("utf-8")
        with open(os.path.join(self.path, name), encoding="utf-8") as f:
            return f.read()

    def gdal_path(self, pol):
        name = self.files[self._check(pol)]["measurement"]
        if self.is_zip:
            return "/vsizip/" + self.path.replace("\\", "/") + "/" + name
        return os.path.join(self.path, name)

    def info(self, pol):
        return parse_product_info(self.read_text(pol, "annotation"))

    def geolocation_grid(self, pol):
        return parse_geolocation_grid(self.read_text(pol, "annotation"))

    def window_for_bbox(self, pol, bbox, margin=200):
        """(col0, row0, width, height) covering a lon/lat box (lon_min, lat_min, lon_max, lat_max)."""
        info = self.info(pol)
        return window_for_bbox(self.geolocation_grid(pol), bbox, info["lines"], info["samples"], margin)

    def read_dn(self, pol, window=None):
        """GRD pixel values (uint16) for a window (col0, row0, width, height), or the whole image."""
        ds = gdal.Open(self.gdal_path(pol))
        if ds is None:
            raise IOError(f"Cannot open {self.gdal_path(pol)}")
        if window is None:
            window = (0, 0, ds.RasterXSize, ds.RasterYSize)
        return ds.GetRasterBand(1).ReadAsArray(*window)

    def backscatter(self, pol, window=None, kind="sigmaNought", remove_noise=True):
        """Calibrated backscatter (linear power, float32) for a window; NaN where there is no data."""
        if window is None:
            info = self.info(pol)
            window = (0, 0, info["samples"], info["lines"])
        dn = self.read_dn(pol, window)
        cal = parse_calibration(self.read_text(pol, "calibration"), kind)
        noise = parse_noise(self.read_text(pol, "noise")) if remove_noise else None
        return calibrate(dn, window[1], window[0], cal, noise)
