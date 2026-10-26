"""Processing provider: the group of tools shown in the Processing Toolbox."""

import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithms.water_map import WaterMapAlgorithm


class SarFloodDetectionProvider(QgsProcessingProvider):
    def loadAlgorithms(self):
        self.addAlgorithm(WaterMapAlgorithm())

    def id(self):
        return "sarflood"

    def name(self):
        return "SAR Flood Detection"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.svg"))
