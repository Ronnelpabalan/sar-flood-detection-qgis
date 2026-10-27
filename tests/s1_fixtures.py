"""Small synthetic Sentinel-1 annotation files for tests (same layout as real products)."""


def _join(values):
    return " ".join(str(v) for v in values)


def calibration_xml(vectors):
    """vectors: list of (line, pixels, sigma0 values); beta0 and gamma get the same values."""
    body = "".join(
        f"<calibrationVector><azimuthTime>2026-10-03T22:50:43</azimuthTime><line>{line}</line>"
        f"<pixel count=\"{len(px)}\">{_join(px)}</pixel>"
        f"<sigmaNought count=\"{len(v)}\">{_join(v)}</sigmaNought>"
        f"<betaNought count=\"{len(v)}\">{_join(v)}</betaNought>"
        f"<gamma count=\"{len(v)}\">{_join(v)}</gamma>"
        f"<dn count=\"{len(v)}\">{_join(v)}</dn></calibrationVector>"
        for line, px, v in vectors)
    return f"<calibration><calibrationVectorList count=\"{len(vectors)}\">{body}</calibrationVectorList></calibration>"


def noise_xml(range_vectors, azimuth_blocks=(), old_format=False):
    """range_vectors: list of (line, pixels, values).
    azimuth_blocks: list of (first_line, last_line, first_sample, last_sample, lines, values).
    old_format: write the pre-IPF-2.9 layout (noiseVector / noiseLut, no azimuth table).
    """
    vec, lut = ("noiseVector", "noiseLut") if old_format else ("noiseRangeVector", "noiseRangeLut")
    rng = "".join(
        f"<{vec}><azimuthTime>2026-10-03T22:50:43</azimuthTime><line>{line}</line>"
        f"<pixel count=\"{len(px)}\">{_join(px)}</pixel><{lut} count=\"{len(v)}\">{_join(v)}</{lut}></{vec}>"
        for line, px, v in range_vectors)
    az = "".join(
        f"<noiseAzimuthVector><swath>IW{i + 1}</swath><firstAzimuthLine>{fl}</firstAzimuthLine>"
        f"<firstRangeSample>{fs}</firstRangeSample><lastAzimuthLine>{ll}</lastAzimuthLine>"
        f"<lastRangeSample>{ls}</lastRangeSample><line count=\"{len(lines)}\">{_join(lines)}</line>"
        f"<noiseAzimuthLut count=\"{len(v)}\">{_join(v)}</noiseAzimuthLut></noiseAzimuthVector>"
        for i, (fl, ll, fs, ls, lines, v) in enumerate(azimuth_blocks))
    if old_format:
        return f"<noise><noiseVectorList count=\"{len(range_vectors)}\">{rng}</noiseVectorList></noise>"
    return (f"<noise><noiseRangeVectorList count=\"{len(range_vectors)}\">{rng}</noiseRangeVectorList>"
            f"<noiseAzimuthVectorList count=\"{len(azimuth_blocks)}\">{az}</noiseAzimuthVectorList></noise>")


def annotation_xml(lines, samples, points, pol="VH"):
    """points: list of (line, pixel, latitude, longitude, height)."""
    grid = "".join(
        f"<geolocationGridPoint><azimuthTime>2026-10-03T22:50:43</azimuthTime>"
        f"<slantRangeTime>0.005</slantRangeTime><line>{ln}</line><pixel>{px}</pixel>"
        f"<latitude>{lat}</latitude><longitude>{lon}</longitude><height>{h}</height>"
        f"<incidenceAngle>39.0</incidenceAngle><elevationAngle>34.0</elevationAngle></geolocationGridPoint>"
        for ln, px, lat, lon, h in points)
    return (
        "<product><adsHeader><missionId>S1D</missionId><productType>GRD</productType>"
        f"<polarisation>{pol}</polarisation><startTime>2026-10-03T22:50:43.258500</startTime>"
        "<stopTime>2026-10-03T22:51:08.257110</stopTime></adsHeader>"
        "<generalAnnotation><productInformation><pass>Descending</pass></productInformation></generalAnnotation>"
        "<imageAnnotation><imageInformation>"
        f"<numberOfLines>{lines}</numberOfLines><numberOfSamples>{samples}</numberOfSamples>"
        "<rangePixelSpacing>1.0e+01</rangePixelSpacing><azimuthPixelSpacing>1.0e+01</azimuthPixelSpacing>"
        "</imageInformation></imageAnnotation>"
        f"<geolocationGrid><geolocationGridPointList count=\"{len(points)}\">{grid}"
        "</geolocationGridPointList></geolocationGrid></product>")


def simple_grid(lines, samples, lon0=105.0, lat0=20.0, step=0.001):
    """Geolocation points on a regular lon/lat grid: pixel -> longitude, line -> latitude (decreasing)."""
    pts = []
    for ln in range(0, lines + 1, max(lines // 4, 1)):
        for px in range(0, samples + 1, max(samples // 4, 1)):
            pts.append((ln, px, lat0 - ln * step, lon0 + px * step, 0.0))
    return pts
