"""
Boundary Unit Operation (Sink)
Exit boundary with fixed conditions
"""

import numpy as np
from .unit_operation import UnitOperation


class Boundary(UnitOperation):
    """
    Boundary sink: absorbs material with fixed P, T
    """
    
    def __init__(self, name, P_fixed, T_fixed, z_fixed, compounds):
        """
        Args:
            name: Unit name
            P_fixed: Fixed pressure (Pa)
            T_fixed: Fixed temperature (K)
            z_fixed: Fixed composition
            compounds: List of component names
        """
        super().__init__(name, n_inlets=1, n_outlets=0)
        
        self.P_fixed = P_fixed
        self.T_fixed = T_fixed
        self.z_fixed = np.array(z_fixed) / np.sum(z_fixed)
        self.compounds = compounds
    
    def calculate(self):
        """Boundary acts as sink, absorbs inlet conditions"""
        self.run_user_code('before')
        
        try:
            inlet = self.inlet_ports[0].stream
            
            if inlet is None:
                self.run_user_code('after')
                return
            
            # Inlet flow is determined by connected unit
            # Boundary just accepts whatever comes in
            # and could impose back-pressure if needed
            
        except Exception as e:
            print(f"Error in {self.name}: {e}")
        
        self.run_user_code('after')