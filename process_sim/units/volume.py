"""
Dynamic Pipe Segment / Holdup Block

Represents the total holdup of all associated piping segments.
Accumulates material from associated streams, maintains internal state,
and allows draining even if inlet flow is zero.
"""

import numpy as np
import sys
from pathlib import Path

# Proper path handling for imports
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from units.unit_operation import UnitOperation

class Volume(UnitOperation):
    """
    Dynamic pipe segment (holdup block).
    Accumulates material from associated streams, maintains internal state.
    """

    def __init__(self, name, V_total, compounds, streams_governed=None, P_init=101325, T_init=298.15, z_init=None):
        """
        Args:
            name: Unit name
            V_total: Total volume (m³)
            compounds: List of component names
            streams_governed: List of Stream objects whose holdup is included in this block
            P_init: Initial pressure (Pa)
            T_init: Initial temperature (K)
            z_init: Initial composition (list or array)
        """
        super().__init__(name, n_inlets=1, n_outlets=1)
        self.V_total = V_total
        self.compounds = compounds
        self.n_compounds = len(compounds)
        self.P = P_init
        self.T = T_init
        self.z = np.array(z_init) / np.sum(z_init) if z_init is not None else np.ones(self.n_compounds) / self.n_compounds

        # Holdup state
        self.m_total = self._init_mass(P_init, T_init, self.z)
        self.streams_governed = streams_governed if streams_governed is not None else []

        # History for plotting
        self.P_history = []
        self.T_history = []
        self.m_history = []
        self.arrayed_vars = {}

    def _init_mass(self, P, T, z):
        """
        Estimate initial mass from ideal gas law.
        """
        try:
            from thermoprops import Chemical
            MW_avg = sum(z[i] * Chemical(self.compounds[i]).MW for i in range(self.n_compounds))
            R = 8.314  # J/mol·K
            n_moles = (P * self.V_total) / (R * 1000 * T)
            return n_moles * MW_avg
        except Exception:
            return 1.0  # fallback

    def calculate(self, dt):
        """
        Dynamic holdup calculation for a pipe segment.
        """
        self.run_user_code('before')
        try:
            inlet = self.inlet_ports[0].stream
            outlet = self.outlet_ports[0].stream

            if inlet is None or outlet is None:
                self.run_user_code('after')
                return

            # --- Mass balance (holdup) ---
            F_in = inlet.F if inlet.F is not None else 0.0
            F_out = outlet.F if outlet.F is not None else F_in

            # If holdup exists, allow draining; if not, outlet flow is zero
            if self.m_total > 0:
                max_possible_out = self.m_total / dt
                F_out = min(F_out, max_possible_out)
            else:
                F_out = 0.0

            # Update holdup
            dm = (F_in - F_out) * dt
            m_prev = self.m_total
            self.m_total = max(self.m_total + dm, 0.0)

            # --- Mixing: update composition and temperature ---
            if self.m_total > 1e-8:
                # Weighted mixing with incoming stream
                self.z = (self.z * m_prev + inlet.z * F_in * dt) / self.m_total
                self.z = self.z / np.sum(self.z)
                self.T = (self.T * m_prev + inlet.T * F_in * dt) / self.m_total
            else:
                # If empty, set to inlet
                self.z = inlet.z.copy()
                self.T = inlet.T

            # --- Pressure calculation (ideal gas) ---
            try:
                from thermoprops import Chemical
                MW_avg = sum(self.z[i] * Chemical(self.compounds[i]).MW for i in range(self.n_compounds))
                n_moles = self.m_total / MW_avg if MW_avg > 0 else 0.0
                R = 8.314
                if self.V_total > 0:
                    self.P = (n_moles * 1000 * R * self.T) / self.V_total
                self.P = max(1e4, min(self.P, 1e8))  # Clamp to 0.1 to 1000 bar
            except Exception:
                pass

            # --- Outlet stream conditions ---
            outlet.F = F_out
            outlet.P = self.P
            outlet.T = self.T
            outlet.z = self.z.copy()

        except Exception as e:
            print(f"Error in {self.name}: {e}")
            import traceback
            traceback.print_exc()
            # Safe fallback
            if outlet is not None and inlet is not None:
                outlet.P = inlet.P
                outlet.T = inlet.T
                outlet.F = 0.0
                outlet.z = inlet.z.copy()
        self.run_user_code('after')

    def store_history(self):
        """
        Store history for plotting.
        """
        self.P_history.append(self.P)
        self.T_history.append(self.T)
        self.m_history.append(self.m_total)
        for var, value in self.arrayed_vars.items():
            hist_name = f"{var}_history"
            if not hasattr(self, hist_name):
                setattr(self, hist_name, [])
            getattr(self, hist_name).append(value)
