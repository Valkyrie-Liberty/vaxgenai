"""
Advanced Modules for VaxGenAI

This package contains advanced modules for vaccine design against incurable diseases.
"""

# Import only the modules that are properly implemented
try:
    from .neoantigen_identification import NeoantigenIdentifier
    __all__ = ['NeoantigenIdentifier']
except ImportError:
    pass

try:
    from .mhc_class_ii_prediction import MHCClassIIPredictor
    __all__.append('MHCClassIIPredictor')
except ImportError:
    pass

try:
    from .immunogenicity_prediction import ImmunogenicityPredictor
    __all__.append('ImmunogenicityPredictor')
except ImportError:
    pass

try:
    from .immune_evasion_modeling import ImmuneEvasionModeler
    __all__.append('ImmuneEvasionModeler')
except ImportError:
    pass

try:
    from .experimental_validation import ExperimentalValidationIntegrator
    __all__.append('ExperimentalValidationIntegrator')
except ImportError:
    pass

try:
    from .conformational_bcell_prediction import ConformationalBCellPredictor
    __all__.append('ConformationalBCellPredictor')
except ImportError:
    pass

# Initialize __all__ if not already defined
if '__all__' not in locals():
    __all__ = []

__version__ = "2.0.0"

