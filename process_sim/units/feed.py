"""
Feed Boundary Unit Operation
Source of material into the system
"""

import numpy as np
from .unit_operation import UnitOperation


class FeedBoundary(UnitOperation):
    """
    Feed boundary: fixed source with P, T, composition, and optional fixed flow
    """
    
    def __init__(self, name, P_fixed, T_fixed, z_fixed, compounds, F_fixed=None):
        """
        Args:
            name: Unit name
            P_fixed: Fixed pressure (Pa)
            T_fixed: Fixed temperature (K)
            z_fixed: Fixed composition array
            compounds: List of component names
            F_fixed: Optional fixed flow rate (kg/s)
        """
        super().__init__(name, n_inlets=0, n_outlets=1)
        
        self.P_fixed = P_fixed
        self.T_fixed = T_fixed
        self.z_fixed = np.array(z_fixed) / np.sum(z_fixed)
        self.compounds = compounds
        self.F_fixed = F_fixed
        self.parameters = {}
    
    def calculate(self):
        """Update outlet with fixed conditions"""
        self.run_user_code('before')
        
        try:
            outlet = self.outlet_ports[0].stream
            
            if outlet is None:
                self.run_user_code('after')
                return
            
            # Set fixed conditions
            outlet.P = self.P_fixed
            outlet.T = self.T_fixed
            outlet.z = self.z_fixed
            
            # If F_fixed is set, use it
            if self.F_fixed is not None:
                outlet.F = self.F_fixed
        
        except Exception as e:
            print(f"Error in {self.name}: {e}")
        
        self.run_user_code('after')