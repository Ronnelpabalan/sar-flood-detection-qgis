"""Plugin entry point: registers the Processing provider with QGIS."""

from qgis.core import QgsApplication

from .provider import SarFloodDetectionProvider


class SarFloodDetectionPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initProcessing(self):
        """Also called by QGIS without a GUI (e.g. qgis_process), because metadata sets hasProcessingProvider."""
        self.provider = SarFloodDetectionProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()

    def unload(self):
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
