"""
Thermodynamic Properties and Flash Calculations Package
Peng-Robinson EOS based calculations using thermo library
"""

# Import from thermo library
from thermo import Mixture, Chemical

# Import from our modules
from .property_extractor import (
    PropertyExtractor,
    create_property_extractor,
    get_all_properties
)

from .flasher import (
    Flasher,
    create_flasher,
    isenthalpic_expansion,
    adiabatic_isentropic_expansion
)

__all__ = [
    # From thermo library
    'Mixture',
    'Chemical',
    
    # From property_extractor
    'PropertyExtractor',
    'create_property_extractor',
    'get_all_properties',
    
    # From flasher
    'Flasher',
    'create_flasher',
    'isenthalpic_expansion',
    'adiabatic_isentropic_expansion',
]

__version__ = '1.0.0'