"""Product information and geolocation grid from a Sentinel-1 annotation XML file."""

import xml.etree.ElementTree as ET

import numpy as np
from scipy.interpolate import griddata


def parse_product_info(xml_text):
    """Basic facts about the image from annotation/s1*.xml."""
    root = ET.fromstring(xml_text)
    return {
        "mission": root.findtext("adsHeader/missionId"),
        "polarisation": root.findtext("adsHeader/polarisation"),
        "start_time": root.findtext("adsHeader/startTime"),
        "stop_time": root.findtext("adsHeader/stopTime"),
        "pass": root.findtext("generalAnnotation/productInformation/pass"),
        "lines": int(root.findtext("imageAnnotation/imageInformation/numberOfLines")),
        "samples": int(root.findtext("imageAnnotation/imageInformation/numberOfSamples")),
        "range_spacing": float(root.findtext("imageAnnotation/imageInformation/rangePixelSpacing")),
        "azimuth_spacing": float(root.findtext("imageAnnotation/imageInformation/azimuthPixelSpacing")),
    }


def parse_geolocation_grid(xml_text):
    """Geolocation grid points (ESA's GCPs) as a dict of arrays.

    Keys: line, pixel, latitude, longitude, height, incidence_angle.
    """
    root = ET.fromstring(xml_text)
    keys = {"line": "line", "pixel": "pixel", "latitude": "latitude", "longitude": "longitude",
            "height": "height", "incidence_angle": "incidenceAngle"}
    points = {k: [] for k in keys}
    for p in root.iter("geolocationGridPoint"):
        for k, tag in keys.items():
            points[k].append(float(p.findtext(tag)))
    if not points["line"]:
        raise ValueError("No geolocationGridPoint elements found.")
    return {k: np.array(v) for k, v in points.items()}


def window_for_bbox(grid, bbox, n_lines, n_samples, margin=200):
    """Image window (col0, row0, width, height) that covers a lon/lat bounding box.

    bbox is (lon_min, lat_min, lon_max, lat_max). Pixel and line positions are
    interpolated linearly from the geolocation grid, which ignores fine terrain
    detail; margin (pixels) adds room for that. Raises ValueError if the box
    does not overlap the image.
    """
    lon_min, lat_min, lon_max, lat_max = bbox
    lon, lat = np.meshgrid(np.linspace(lon_min, lon_max, 21), np.linspace(lat_min, lat_max, 21))
    known = np.column_stack([grid["longitude"], grid["latitude"]])
    query = np.column_stack([lon.ravel(), lat.ravel()])
    pix = griddata(known, grid["pixel"], query, method="linear")
    lin = griddata(known, grid["line"], query, method="linear")
    inside = ((grid["longitude"] >= lon_min) & (grid["longitude"] <= lon_max)
              & (grid["latitude"] >= lat_min) & (grid["latitude"] <= lat_max))
    pix = np.concatenate([pix[np.isfinite(pix)], grid["pixel"][inside]])
    lin = np.concatenate([lin[np.isfinite(lin)], grid["line"][inside]])
    if pix.size == 0:
        raise ValueError("The area does not overlap the image.")

    eps = 1e-6                                       # so 99.9999999 counts as 100
    col0 = max(int(np.floor(pix.min() + eps)) - margin, 0)
    row0 = max(int(np.floor(lin.min() + eps)) - margin, 0)
    col1 = min(int(np.ceil(pix.max() - eps)) + margin, n_samples)
    row1 = min(int(np.ceil(lin.max() - eps)) + margin, n_lines)
    if col1 <= col0 or row1 <= row0:
        raise ValueError("The area does not overlap the image.")
    return col0, row0, col1 - col0, row1 - row0
