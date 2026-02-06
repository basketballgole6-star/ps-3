"""
Advanced Valve Unit Operation
Includes proper thermodynamic calculations with:
- Isenthalpic expansion (throttling)
- Compressibility correction (Z-factor)
- Two-phase flow handling with flash calculations
- Dynamic performance curves
- Real boundary pressure calculations
"""

import numpy as np
import sys
import os
from pathlib import Path

# Proper path handling for imports
current_dir = Path(__file__).parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

# Import from same package
from units.unit_operation import UnitOperation
from thermoprops import PropertyExtractor, Flasher, isenthalpic_expansion


class Valve(UnitOperation):
    """
    Advanced Valve with:
    - Isenthalpic expansion (throttling)
    - Compressibility correction (Z-factor)
    - Two-phase flow handling with flash calculations
    - Dynamic performance curves
    - Real boundary pressure calculations
    """
    
    def __init__(self, name, Cv=10.0, opening=1.0, boundary_P=None, valve_type='globe'):
        """
        Args:
            name (str): Valve name
            Cv (float): Valve flow coefficient (US gallons/min per sqrt(psi))
            opening (float): Valve opening fraction (0-1)
            boundary_P (float): Fixed downstream boundary pressure (Pa)
            valve_type (str): 'globe', 'ball', 'butterfly', 'gate'
        """
        super().__init__(name, n_inlets=1, n_outlets=1)
        
        self.valve_type = valve_type
        self.parameters = {
            'Cv': Cv,
            'opening': opening,
            'boundary_P': boundary_P,
        }
        
        self.performance_curves = self._get_performance_curve()
        self.results = {}
    
    def _get_performance_curve(self):
        """Get valve performance curve based on valve type"""
        if self.valve_type == 'globe':
            return lambda opening: opening ** 2
        elif self.valve_type == 'ball':
            return lambda opening: 0.1 * opening + 0.9 * opening ** 2
        elif self.valve_type == 'butterfly':
            return lambda opening: opening ** 1.5 if opening > 0.3 else 0.0
        elif self.valve_type == 'gate':
            return lambda opening: opening if opening > 0.1 else 0.0
        else:
            return lambda opening: opening
    
    def _get_inlet_properties(self):
        """Get inlet stream properties"""
        inlet_stream = self.inlet_ports[0].stream
        if inlet_stream is None:
            return None
        
        # Convert Stream.z array to dict for thermoprops
        composition_dict = {
            inlet_stream.compounds[i]: inlet_stream.z[i] 
            for i in range(len(inlet_stream.compounds))
        }
        
        return {
            'P': inlet_stream.P,
            'T': inlet_stream.T,
            'F': inlet_stream.F,
            'components': inlet_stream.compounds,
            'composition': composition_dict,
        }
    
    def _calculate_outlet_conditions_isenthalpic(self, P_outlet):
        """Calculate outlet conditions using isenthalpic expansion"""
        try:
            inlet_props = self._get_inlet_properties()
            if inlet_props is None:
                return None
            
            flash_result = isenthalpic_expansion(
                components=inlet_props['components'],
                composition=inlet_props['composition'],
                inlet_pressure=inlet_props['P'],
                inlet_temperature=inlet_props['T'],
                outlet_pressure=P_outlet
            )
            
            if flash_result is None:
                return None
            
            self.results['flash_result'] = flash_result
            self.results['outlet_temp'] = flash_result['temperature']
            self.results['vapor_fraction'] = flash_result['vapor_fraction']
            self.results['outlet_phase'] = flash_result['phase']
            
            return flash_result
        
        except Exception as e:
            print(f"Error in isenthalpic expansion: {e}")
            return None
    
    def _calculate_boundary_pressure(self, P_inlet, dP_target):
        """Calculate boundary pressure with compressibility correction"""
        try:
            inlet_props = self._get_inlet_properties()
            if inlet_props is None:
                return P_inlet - dP_target
            
            inlet_extractor = PropertyExtractor(
                inlet_props['components'],
                inlet_props['composition'],
                inlet_props['P'],
                inlet_props['T']
            )
            Z_inlet = inlet_extractor.get_compressibility_factor() or 1.0
            
            P_outlet = P_inlet - dP_target
            
            for iteration in range(3):
                outlet_extractor = PropertyExtractor(
                    inlet_props['components'],
                    inlet_props['composition'],
                    P_outlet,
                    inlet_props['T']
                )
                Z_outlet = outlet_extractor.get_compressibility_factor() or 1.0
                
                Z_avg = (Z_inlet + Z_outlet) / 2
                dP_corrected = dP_target * (1.0 / Z_avg if Z_avg > 0 else 1.0)
                
                P_outlet = P_inlet - dP_corrected
                
                if P_outlet < 1e5:
                    P_outlet = 1e5
                    break
            
            self.results['Z_inlet'] = Z_inlet
            self.results['Z_outlet'] = Z_outlet
            self.results['P_outlet_calculated'] = P_outlet
            
            return P_outlet
        
        except Exception as e:
            print(f"Error calculating boundary pressure: {e}")
            return P_inlet - dP_target

    def _calculate_compressible_flow(self, Cv_effective, dP):
        """Calculate compressible flow - FIXED VERSION"""
        try:
            inlet_props = self._get_inlet_properties()
            if inlet_props is None:
                return 0.0
            
            # Get properties
            P_inlet = inlet_props['P']
            P_outlet = P_inlet - dP
            
            # Get molecular weight
            inlet_extractor = PropertyExtractor(
                inlet_props['components'],
                inlet_props['composition'],
                P_inlet,
                inlet_props['T']
            )
            MW = inlet_extractor.get_molar_weight()  # kg/kmol
            if MW is None or MW == 0:
                return 0.0
            
            # Get outlet conditions
            outlet_props = self._calculate_outlet_conditions_isenthalpic(P_outlet)
            if outlet_props is None:
                return 0.0
            
            T_outlet = outlet_props['temperature']
            Z_outlet = inlet_extractor.get_compressibility_factor() or 1.0
            
            # Convert Cv to SI (m³/s at 1 bar differential with water)
            # Cv in US gallons/min/√psi → multiply by 0.0000241 to get m³/s/√Pa
            Cv_SI = Cv_effective * 0.0000241
            
            # Compressible flow equation (gas)
            # W = Cv * Y * sqrt(ρ_inlet * ΔP)
            # where Y = expansion factor ≈ 1 - (ΔP)/(3*P_inlet) for subsonic
            
            # Get inlet density
            rho_inlet = inlet_extractor.get_density()  # kg/m³
            if rho_inlet is None or rho_inlet <= 0:
                # Fallback: ideal gas density
                R_specific = 8.314 / MW * 1000  # J/(kg·K)
                rho_inlet = P_inlet / (R_specific * inlet_props['T'])
            
            # Expansion factor (simplified)
            Y = 1.0 - dP / (3.0 * P_inlet)
            Y = max(0.5, min(Y, 1.0))  # Limit to reasonable range
            
            # Mass flow rate (kg/s)
            # F = Cv * Y * sqrt(ρ * ΔP)
            F = Cv_SI * Y * np.sqrt(rho_inlet * dP)
            
            # Store results
            self.results['mass_flow'] = F
            self.results['Cv_effective'] = Cv_effective
            self.results['dP'] = dP
            self.results['Y_expansion'] = Y
            self.results['rho_inlet'] = rho_inlet
            
            return F
            
        except Exception as e:
            print(f"Error calculating compressible flow: {e}")
            return 0.0

    def _check_choked_flow(self, P_inlet, P_outlet):
        """Check if flow is choked"""
        try:
            inlet_props = self._get_inlet_properties()
            if inlet_props is None:
                return False
            
            inlet_extractor = PropertyExtractor(
                inlet_props['components'],
                inlet_props['composition'],
                P_inlet,
                inlet_props['T']
            )
            gamma = inlet_extractor.get_gamma()
            
            if gamma is None or gamma <= 1.0:
                return False
            
            critical_pressure_ratio = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
            actual_ratio = P_outlet / P_inlet if P_inlet > 0 else 1.0
            
            is_choked = actual_ratio <= critical_pressure_ratio
            
            self.results['gamma'] = gamma
            self.results['critical_pressure_ratio'] = critical_pressure_ratio
            self.results['actual_pressure_ratio'] = actual_ratio
            self.results['is_choked'] = is_choked
            
            return is_choked
        
        except Exception as e:
            print(f"Error checking choked flow: {e}")
            return False

    def calculate(self):
        """Main valve calculation"""
        self.run_user_code('before')

        try:
            inlet = self.inlet_ports[0].stream
            outlet = self.outlet_ports[0].stream
            
            if inlet is None or outlet is None:
                self.results['status'] = 'error'
                self.run_user_code('after')
                return
            
            opening = self.parameters['opening']
            if opening < 0.0001:
                # FIXED: Set all outlet properties when valve is closed
                outlet.F = 0.0
                outlet.P = inlet.P
                outlet.T = inlet.T
                outlet.z = inlet.z.copy()
                inlet.F = 0.0
                self.results['status'] = 'closed'
                self.run_user_code('after')
                return
            
            inlet_props = self._get_inlet_properties()
            if inlet_props is None:
                # FIXED: Set fallback values
                outlet.F = 0.0
                outlet.P = inlet.P
                outlet.T = inlet.T
                outlet.z = inlet.z.copy()
                inlet.F = 0.0
                self.results['status'] = 'error_no_inlet_props'
                self.run_user_code('after')
                return
            
            Cv_base = self.parameters['Cv']
            boundary_P = self.parameters['boundary_P']
            
            Cv_multiplier = self.performance_curves(opening)
            Cv_effective = Cv_base * Cv_multiplier
            
            if boundary_P is not None:
                P_outlet = boundary_P
                self.results['boundary_pressure_type'] = 'fixed'
            else:
                dP_estimated = 0.3 * inlet_props['P']
                P_outlet = self._calculate_boundary_pressure(inlet_props['P'], dP_estimated)
                self.results['boundary_pressure_type'] = 'calculated'
            
            dP = inlet_props['P'] - P_outlet
            
            is_choked = self._check_choked_flow(inlet_props['P'], P_outlet)
            
            if is_choked:
                self.results['warning'] = 'choked_flow'
            
            if dP > 1e3:
                F = self._calculate_compressible_flow(Cv_effective, dP)
            else:
                F = 0.0
            
            # FIXED: Handle None return from flash
            outlet_conditions = self._calculate_outlet_conditions_isenthalpic(P_outlet)
            
            if outlet_conditions is None:
                # FIXED: Set fallback values - isenthalpic means no temp change as fallback
                outlet.F = F if F is not None else 0.0
                outlet.P = P_outlet
                outlet.T = inlet.T  # Fallback: assume no temperature change
                outlet.z = np.array([inlet_props['composition'].get(comp, 0) for comp in inlet.compounds])
                inlet.F = F if F is not None else 0.0
                self.results['status'] = 'error_flash_failed'
                self.results['calculation_valid'] = False
                self.run_user_code('after')
                return
            
            # FIXED: Check that temperature is valid
            outlet_temp = outlet_conditions.get('temperature')
            if outlet_temp is None or outlet_temp <= 0:
                outlet_temp = inlet.T  # Fallback
            
            # Set outlet conditions
            outlet.F = F if F is not None else 0.0
            outlet.P = P_outlet
            outlet.T = outlet_temp
            outlet.z = np.array([inlet_props['composition'].get(comp, 0) for comp in inlet.compounds])
            
            # Set inlet flow (valve determines flow)
            inlet.F = F if F is not None else 0.0
            
            self.results['status'] = 'open'
            self.results['calculation_valid'] = True
            
        except Exception as e:
            print(f"Error in valve {self.name} calculation: {e}")
            import traceback
            traceback.print_exc()
            # FIXED: Set safe fallback values even on exception
            if inlet is not None and outlet is not None:
                outlet.F = 0.0
                outlet.P = inlet.P
                outlet.T = inlet.T
                outlet.z = inlet.z.copy()
                inlet.F = 0.0
            self.results['status'] = 'error'
            self.results['calculation_valid'] = False
        
        self.run_user_code('after')

    def get_results(self):
        """Return calculation results"""
        return self.results