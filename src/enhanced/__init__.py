"""
Enhanced modules for VaxGenAI

This package contains enhanced implementations of the core VaxGenAI modules
to address the limitations of the original prototype.
"""

__version__ = "0.2.0"

from .epitope_prediction import (
    EnhancedEpitopePrediction,
    AlphaFoldIntegration,
    PopulationCoverageAnalysis
)
from .safety_filter import EnhancedSafetyFilter

