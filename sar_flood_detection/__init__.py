"""SAR Flood Detection: QGIS Processing plugin for Sentinel-1 flood mapping."""


def classFactory(iface):
    """Called by QGIS to create the plugin."""
    from .plugin import SarFloodDetectionPlugin
    return SarFloodDetectionPlugin(iface)
