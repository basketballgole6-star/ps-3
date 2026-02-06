"""
Equation System Formulation for Process Simulation
Builds and manages all equations (mass, energy, momentum, constitutive)
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Callable
import logging

logger = logging.getLogger(__name__)


class EquationSystem:
    """
    Formulates the complete equation system for process simulation.
    
    Equation types:
    1. Differential equations: dm/dt, dU/dt (volumes)
    2. Algebraic equations: valve characteristics, stream properties, boundaries
    3. Constraint equations: flow continuity, phase equilibrium
    """
    
    def __init__(self, simulator):
        """
        Args:
            simulator: DynamicSimulator instance containing all units and streams
        """
        self.simulator = simulator
        
        # Equation components
        self.differential_equations = []
        self.algebraic_equations = []
        self.constraint_equations = []
        
        # Variable mapping
        self.differential_vars = []  # List of (unit, var_name)
        self.algebraic_vars = []     # List of (stream/unit, var_name)
        
        # Jacobian structure
        self.jacobian_sparsity = None
        
        # Cache
        self._equations_built = False
        
    def build_equations(self):
        """
        Build complete equation system from simulator units and streams
        """
        logger.info("Building equation system...")
        
        # Clear existing
        self.differential_equations = []
        self.algebraic_equations = []
        self.constraint_equations = []
        self.differential_vars = []
        self.algebraic_vars = []
        
        # 1. Identify variables
        self._identify_variables()
        
        # 2. Build differential equations (volumes)
        self._build_volume_equations()
        
        # 3. Build algebraic equations (valves, streams)
        self._build_valve_equations()
        self._build_stream_equations()
        self._build_boundary_equations()
        
        # 4. Build constraint equations
        self._build_flow_continuity_equations()
        
        self._equations_built = True
        
        n_diff = len(self.differential_equations)
        n_alg = len(self.algebraic_equations)
        n_const = len(self.constraint_equations)
        
        logger.info(f"Equation system built: {n_diff} differential, "
                   f"{n_alg} algebraic, {n_const} constraints")
        logger.info(f"Variables: {len(self.differential_vars)} differential, "
                   f"{len(self.algebraic_vars)} algebraic")
        
    def _identify_variables(self):
        """
        Identify all variables in the system
        
        Differential variables (state variables):
        - Volume mass (m_total)
        - Volume internal energy (U_total)
        - Volume component moles (n_i)
        
        Algebraic variables:
        - Stream pressures (P)
        - Stream temperatures (T)
        - Stream flows (F)
        - Volume pressures (P)
        - Volume temperatures (T)
        """
        # Differential variables from volumes
        for unit in self.simulator.units:
            if hasattr(unit, 'P_history'):  # Is a volume
                self.differential_vars.append((unit, 'm_total'))
                self.differential_vars.append((unit, 'U_total'))
                
                # Component moles
                if hasattr(unit, 'n_moles'):
                    for i in range(len(unit.n_moles)):
                        self.differential_vars.append((unit, f'n_moles_{i}'))
        
        # Algebraic variables from streams
        for stream_name, stream in self.simulator.streams.items():
            self.algebraic_vars.append((stream, 'P'))
            self.algebraic_vars.append((stream, 'T'))
            self.algebraic_vars.append((stream, 'F'))
        
        # Algebraic variables from volumes (P, T are algebraic, computed from state)
        for unit in self.simulator.units:
            if hasattr(unit, 'P_history'):
                self.algebraic_vars.append((unit, 'P'))
                self.algebraic_vars.append((unit, 'T'))
    
    def _build_volume_equations(self):
        """
        Build differential equations for volumes
        
        For each volume:
        1. dm/dt = F_in - F_out
        2. dU/dt = F_in*H_in - F_out*H - Q
        3. dn_i/dt = F_in*z_in[i]/MW_in - F_out*z[i]/MW_out
        """
        for unit in self.simulator.units:
            if not hasattr(unit, 'P_history'):  # Not a volume
                continue
            
            # Mass balance
            def mass_balance(t, state_vars, alg_vars):
                """dm/dt = F_in - F_out"""
                F_in = self._get_inlet_flow(unit, alg_vars)
                F_out = self._get_outlet_flow(unit, alg_vars)
                return F_in - F_out
            
            self.differential_equations.append({
                'name': f'{unit.name}_mass_balance',
                'unit': unit,
                'variable': 'm_total',
                'function': mass_balance,
                'type': 'mass_balance'
            })
            
            # Energy balance
            def energy_balance(t, state_vars, alg_vars):
                """dU/dt = F_in*H_in - F_out*H - Q"""
                F_in = self._get_inlet_flow(unit, alg_vars)
                F_out = self._get_outlet_flow(unit, alg_vars)
                H_in = self._get_inlet_enthalpy(unit, alg_vars)
                H_out = self._get_outlet_enthalpy(unit, alg_vars)
                
                # Heat transfer
                Q = 0.0
                if hasattr(unit, 'UA'):
                    T_volume = self._get_volume_temp(unit, alg_vars)
                    Q = unit.UA * (T_volume - unit.T_ambient)
                
                return F_in * H_in - F_out * H_out - Q
            
            self.differential_equations.append({
                'name': f'{unit.name}_energy_balance',
                'unit': unit,
                'variable': 'U_total',
                'function': energy_balance,
                'type': 'energy_balance'
            })
            
            # Component mass balances
            if hasattr(unit, 'n_moles'):
                for i in range(len(unit.n_moles)):
                    def component_balance(t, state_vars, alg_vars, component_idx=i):
                        """dn_i/dt = n_dot_in*z_in[i] - n_dot_out*z[i]"""
                        F_in = self._get_inlet_flow(unit, alg_vars)
                        F_out = self._get_outlet_flow(unit, alg_vars)
                        
                        z_in = self._get_inlet_composition(unit, alg_vars)
                        z_out = self._get_outlet_composition(unit, alg_vars)
                        
                        MW_in = self._get_mixture_MW(z_in, unit.compounds)
                        MW_out = self._get_mixture_MW(z_out, unit.compounds)
                        
                        n_dot_in = F_in / MW_in * 1000 if MW_in > 0 else 0  # kmol/s
                        n_dot_out = F_out / MW_out * 1000 if MW_out > 0 else 0
                        
                        return n_dot_in * z_in[component_idx] - n_dot_out * z_out[component_idx]
                    
                    self.differential_equations.append({
                        'name': f'{unit.name}_component_{i}_balance',
                        'unit': unit,
                        'variable': f'n_moles_{i}',
                        'function': component_balance,
                        'type': 'component_balance'
                    })
    
    def _build_valve_equations(self):
        """
        Build algebraic equations for valves
        
        For each valve:
        1. Flow equation: F = Cv * Y * sqrt(rho * dP)
        2. Outlet enthalpy: H_out = H_in (isenthalpic)
        3. Outlet composition: z_out = z_in
        """
        from units.valve import Valve
        
        for unit in self.simulator.units:
            if not isinstance(unit, Valve):
                continue
            
            # Valve characteristic equation
            def valve_flow_equation(t, state_vars, alg_vars):
                """F_actual - F_calculated = 0"""
                inlet_stream = unit.inlet_ports[0].stream
                outlet_stream = unit.outlet_ports[0].stream
                
                if inlet_stream is None or outlet_stream is None:
                    return 0.0
                
                # Get variables from alg_vars
                P_in = self._get_stream_var(inlet_stream, 'P', alg_vars)
                P_out = self._get_stream_var(outlet_stream, 'P', alg_vars)
                F_actual = self._get_stream_var(outlet_stream, 'F', alg_vars)
                
                # Calculate flow from pressure drop
                dP = P_in - P_out
                if dP < 0:
                    dP = 0
                
                # Valve equation (simplified)
                opening = unit.parameters.get('opening', 1.0)
                Cv = unit.parameters.get('Cv', 10.0)
                
                if opening < 0.0001:
                    F_calculated = 0.0
                else:
                    # Get density
                    T_in = self._get_stream_var(inlet_stream, 'T', alg_vars)
                    rho = self._calculate_density(inlet_stream, P_in, T_in)
                    
                    # Expansion factor
                    Y = 1.0 - dP / (3.0 * P_in) if P_in > 0 else 1.0
                    Y = max(0.5, min(Y, 1.0))
                    
                    # Flow equation
                    Cv_SI = Cv * 0.0000241  # Convert to SI
                    F_calculated = Cv_SI * opening**2 * Y * np.sqrt(rho * dP)
                
                return F_actual - F_calculated
            
            self.algebraic_equations.append({
                'name': f'{unit.name}_flow_equation',
                'unit': unit,
                'function': valve_flow_equation,
                'type': 'valve_characteristic'
            })
            
            # Isenthalpic constraint
            def isenthalpic_constraint(t, state_vars, alg_vars):
                """T_out - T_flash(P_out, H_in) = 0"""
                inlet_stream = unit.inlet_ports[0].stream
                outlet_stream = unit.outlet_ports[0].stream
                
                if inlet_stream is None or outlet_stream is None:
                    return 0.0
                
                # This would require flash calculation
                # For now, simplified: small temperature drop
                T_in = self._get_stream_var(inlet_stream, 'T', alg_vars)
                T_out = self._get_stream_var(outlet_stream, 'T', alg_vars)
                P_in = self._get_stream_var(inlet_stream, 'P', alg_vars)
                P_out = self._get_stream_var(outlet_stream, 'P', alg_vars)
                
                # Simplified: dT ~ 0.1 * dP/P
                dP = P_in - P_out
                dT_expected = 0.1 * dP / P_in * T_in if P_in > 0 else 0
                T_expected = T_in - dT_expected
                
                return T_out - T_expected
            
            self.algebraic_equations.append({
                'name': f'{unit.name}_temperature_constraint',
                'unit': unit,
                'function': isenthalpic_constraint,
                'type': 'thermodynamic_constraint'
            })
    
    def _build_stream_equations(self):
        """
        Build algebraic equations for streams
        
        For each stream:
        1. Composition sum: sum(z_i) = 1
        2. Property consistency
        """
        for stream_name, stream in self.simulator.streams.items():
            
            # Composition constraint
            def composition_constraint(t, state_vars, alg_vars, s=stream):
                """sum(z_i) - 1 = 0"""
                return np.sum(s.z) - 1.0
            
            self.algebraic_equations.append({
                'name': f'{stream_name}_composition_sum',
                'stream': stream,
                'function': composition_constraint,
                'type': 'composition_constraint'
            })
    
    def _build_boundary_equations(self):
        """
        Build equations for boundary conditions (Feed, Product)
        """
        from units.feed import FeedBoundary
        from units.boundary import Boundary
        
        for unit in self.simulator.units:
            # Feed boundary
            if isinstance(unit, FeedBoundary):
                outlet_stream = unit.outlet_ports[0].stream
                
                # Fixed pressure
                def feed_pressure(t, state_vars, alg_vars, u=unit, s=outlet_stream):
                    P = self._get_stream_var(s, 'P', alg_vars)
                    return P - u.P_fixed
                
                self.algebraic_equations.append({
                    'name': f'{unit.name}_pressure_boundary',
                    'unit': unit,
                    'function': feed_pressure,
                    'type': 'boundary_condition'
                })
                
                # Fixed temperature
                def feed_temperature(t, state_vars, alg_vars, u=unit, s=outlet_stream):
                    T = self._get_stream_var(s, 'T', alg_vars)
                    return T - u.T_fixed
                
                self.algebraic_equations.append({
                    'name': f'{unit.name}_temperature_boundary',
                    'unit': unit,
                    'function': feed_temperature,
                    'type': 'boundary_condition'
                })
                
                # Fixed flow
                if hasattr(unit, 'F_fixed') and unit.F_fixed is not None:
                    def feed_flow(t, state_vars, alg_vars, u=unit, s=outlet_stream):
                        F = self._get_stream_var(s, 'F', alg_vars)
                        return F - u.F_fixed
                    
                    self.algebraic_equations.append({
                        'name': f'{unit.name}_flow_boundary',
                        'unit': unit,
                        'function': feed_flow,
                        'type': 'boundary_condition'
                    })
            
            # Product boundary
            elif isinstance(unit, Boundary):
                inlet_stream = unit.inlet_ports[0].stream
                
                # Fixed pressure
                def product_pressure(t, state_vars, alg_vars, u=unit, s=inlet_stream):
                    P = self._get_stream_var(s, 'P', alg_vars)
                    return P - u.P_fixed
                
                self.algebraic_equations.append({
                    'name': f'{unit.name}_pressure_boundary',
                    'unit': unit,
                    'function': product_pressure,
                    'type': 'boundary_condition'
                })
    
    def _build_flow_continuity_equations(self):
        """
        Build flow continuity equations
        
        At each connection point: F_in = F_out (unless accumulating)
        """
        # For non-volume units, flow is continuous
        for unit in self.simulator.units:
            if hasattr(unit, 'P_history'):  # Is a volume - has accumulation
                continue
            
            # Skip boundaries
            from units.feed import FeedBoundary
            from units.boundary import Boundary
            if isinstance(unit, (FeedBoundary, Boundary)):
                continue
            
            # For valves, etc: F_in = F_out
            if len(unit.inlet_ports) > 0 and len(unit.outlet_ports) > 0:
                def flow_continuity(t, state_vars, alg_vars, u=unit):
                    F_in = sum(self._get_stream_var(port.stream, 'F', alg_vars) 
                              for port in u.inlet_ports if port.stream)
                    F_out = sum(self._get_stream_var(port.stream, 'F', alg_vars)
                               for port in u.outlet_ports if port.stream)
                    return F_in - F_out
                
                self.constraint_equations.append({
                    'name': f'{unit.name}_flow_continuity',
                    'unit': unit,
                    'function': flow_continuity,
                    'type': 'flow_continuity'
                })
    
    # Helper methods
    def _get_inlet_flow(self, unit, alg_vars) -> float:
        """Get total inlet flow to unit"""
        F_total = 0.0
        for port in unit.inlet_ports:
            if port.stream:
                F_total += self._get_stream_var(port.stream, 'F', alg_vars)
        return F_total
    
    def _get_outlet_flow(self, unit, alg_vars) -> float:
        """Get total outlet flow from unit"""
        F_total = 0.0
        for port in unit.outlet_ports:
            if port.stream:
                F_total += self._get_stream_var(port.stream, 'F', alg_vars)
        return F_total
    
    def _get_inlet_enthalpy(self, unit, alg_vars) -> float:
        """Get inlet specific enthalpy (J/kg)"""
        if len(unit.inlet_ports) == 0 or not unit.inlet_ports[0].stream:
            return 0.0
        
        stream = unit.inlet_ports[0].stream
        P = self._get_stream_var(stream, 'P', alg_vars)
        T = self._get_stream_var(stream, 'T', alg_vars)
        
        # Calculate enthalpy
        from thermoprops import PropertyExtractor
        comp_dict = {stream.compounds[i]: stream.z[i] for i in range(len(stream.compounds))}
        props = PropertyExtractor(stream.compounds, comp_dict, P, T)
        H = props.get_enthalpy_mass()  # J/kg
        return H if H is not None else 0.0
    
    def _get_outlet_enthalpy(self, unit, alg_vars) -> float:
        """Get outlet specific enthalpy (J/kg)"""
        if len(unit.outlet_ports) == 0 or not unit.outlet_ports[0].stream:
            return 0.0
        
        # Use volume properties if available
        if hasattr(unit, 'P_history'):
            P = self._get_volume_var(unit, 'P', alg_vars)
            T = self._get_volume_var(unit, 'T', alg_vars)
            z = unit.z
        else:
            stream = unit.outlet_ports[0].stream
            P = self._get_stream_var(stream, 'P', alg_vars)
            T = self._get_stream_var(stream, 'T', alg_vars)
            z = stream.z
        
        from thermoprops import PropertyExtractor
        comp_dict = {unit.compounds[i]: z[i] for i in range(len(unit.compounds))}
        props = PropertyExtractor(unit.compounds, comp_dict, P, T)
        H = props.get_enthalpy_mass()
        return H if H is not None else 0.0
    
    def _get_inlet_composition(self, unit, alg_vars) -> np.ndarray:
        """Get inlet composition"""
        if len(unit.inlet_ports) == 0 or not unit.inlet_ports[0].stream:
            return np.ones(len(unit.compounds)) / len(unit.compounds)
        return unit.inlet_ports[0].stream.z
    
    def _get_outlet_composition(self, unit, alg_vars) -> np.ndarray:
        """Get outlet composition"""
        if hasattr(unit, 'z'):
            return unit.z
        if len(unit.outlet_ports) == 0 or not unit.outlet_ports[0].stream:
            return np.ones(len(unit.compounds)) / len(unit.compounds)
        return unit.outlet_ports[0].stream.z
    
    def _get_volume_temp(self, unit, alg_vars) -> float:
        """Get volume temperature"""
        return self._get_volume_var(unit, 'T', alg_vars)
    
    def _get_stream_var(self, stream, var_name: str, alg_vars) -> float:
        """Get stream variable from algebraic variables vector"""
        # Find index in algebraic variables
        for i, (obj, vname) in enumerate(self.algebraic_vars):
            if obj is stream and vname == var_name:
                if alg_vars is not None and i < len(alg_vars):
                    return alg_vars[i]
                else:
                    # Use current value
                    return getattr(stream, var_name)
        
        # Fallback
        return getattr(stream, var_name, 0.0)
    
    def _get_volume_var(self, unit, var_name: str, alg_vars) -> float:
        """Get volume variable from algebraic variables vector"""
        for i, (obj, vname) in enumerate(self.algebraic_vars):
            if obj is unit and vname == var_name:
                if alg_vars is not None and i < len(alg_vars):
                    return alg_vars[i]
                else:
                    return getattr(unit, var_name)
        
        return getattr(unit, var_name, 0.0)
    
    def _get_mixture_MW(self, z: np.ndarray, compounds: List[str]) -> float:
        """Calculate mixture molecular weight"""
        from thermoprops import Chemical
        MW = sum(z[i] * Chemical(compounds[i]).MW for i in range(len(compounds)))
        return MW if MW > 0 else 1.0
    
    def _calculate_density(self, stream, P: float, T: float) -> float:
        """Calculate stream density"""
        from thermoprops import PropertyExtractor
        comp_dict = {stream.compounds[i]: stream.z[i] for i in range(len(stream.compounds))}
        props = PropertyExtractor(stream.compounds, comp_dict, P, T)
        rho = props.get_density()
        return rho if rho is not None and rho > 0 else 1.0
    
    # Evaluation methods
    def get_state_vector(self) -> np.ndarray:
        """Get current state vector (differential variables)"""
        state = []
        for unit, var_name in self.differential_vars:
            if '_' in var_name and var_name.startswith('n_moles'):
                # Component mole
                idx = int(var_name.split('_')[-1])
                state.append(unit.n_moles[idx])
            else:
                state.append(getattr(unit, var_name, 0.0))
        return np.array(state)
    
    def set_state_vector(self, state: np.ndarray):
        """Set state vector (differential variables)"""
        for i, (unit, var_name) in enumerate(self.differential_vars):
            if '_' in var_name and var_name.startswith('n_moles'):
                idx = int(var_name.split('_')[-1])
                unit.n_moles[idx] = state[i]
            else:
                setattr(unit, var_name, state[i])
    
    def get_algebraic_vector(self) -> np.ndarray:
        """Get current algebraic variable vector"""
        alg = []
        for obj, var_name in self.algebraic_vars:
            alg.append(getattr(obj, var_name, 0.0))
        return np.array(alg)
    
    def set_algebraic_vector(self, alg: np.ndarray):
        """Set algebraic variable vector"""
        for i, (obj, var_name) in enumerate(self.algebraic_vars):
            setattr(obj, var_name, alg[i])
    
    def evaluate_differential_equations(self, t: float, state_vars: np.ndarray, 
                                       alg_vars: np.ndarray) -> np.ndarray:
        """
        Evaluate all differential equations: dy/dt = f(t, y, z)
        
        Returns:
            Array of derivatives
        """
        derivatives = np.zeros(len(self.differential_vars))
        
        for i, eq in enumerate(self.differential_equations):
            try:
                derivatives[i] = eq['function'](t, state_vars, alg_vars)
            except Exception as e:
                logger.error(f"Error evaluating {eq['name']}: {e}")
                derivatives[i] = 0.0
        
        return derivatives
    
    def evaluate_algebraic_equations(self, t: float, state_vars: np.ndarray,
                                     alg_vars: np.ndarray) -> np.ndarray:
        """
        Evaluate all algebraic equations: 0 = g(t, y, z)
        
        Returns:
            Array of residuals
        """
        n_alg = len(self.algebraic_equations)
        n_const = len(self.constraint_equations)
        residuals = np.zeros(n_alg + n_const)
        
        # Algebraic equations
        for i, eq in enumerate(self.algebraic_equations):
            try:
                residuals[i] = eq['function'](t, state_vars, alg_vars)
            except Exception as e:
                logger.error(f"Error evaluating {eq['name']}: {e}")
                residuals[i] = 0.0
        
        # Constraint equations
        for i, eq in enumerate(self.constraint_equations):
            try:
                residuals[n_alg + i] = eq['function'](t, state_vars, alg_vars)
            except Exception as e:
                logger.error(f"Error evaluating {eq['name']}: {e}")
                residuals[n_alg + i] = 0.0
        
        return residuals
    
    def get_equation_names(self) -> List[str]:
        """Get names of all equations"""
        names = []
        names.extend([eq['name'] for eq in self.differential_equations])
        names.extend([eq['name'] for eq in self.algebraic_equations])
        names.extend([eq['name'] for eq in self.constraint_equations])
        return names
    
    def get_variable_names(self) -> Tuple[List[str], List[str]]:
        """Get names of all variables"""
        diff_names = [f"{unit.name}.{var}" for unit, var in self.differential_vars]
        alg_names = [f"{obj.name}.{var}" for obj, var in self.algebraic_vars]
        return diff_names, alg_names


if __name__ == "__main__":
    # Test with simple system
    print("EquationSystem module loaded successfully")
    print("Use with DynamicSimulator to build equation system")