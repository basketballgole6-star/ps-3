# units/stream.py
from thermo import Chemical, Mixture
import config.sim_settings as sim_settings

import numpy as np

class Stream:
    """
    Represents a material stream connecting unit operation ports.
    Stores pressure, temperature, flow, and composition.
    """
    def __init__(self, name, compounds, P_init=101325, T_init=298.15, z_init=None, F_init=0.0):
        self.name = name
        self.compounds = compounds
        self.n_compounds = len(compounds)
        
        self.P = P_init
        self.T = T_init
        self.F = F_init
        self.z = np.array(z_init) / np.sum(z_init) if z_init is not None else np.ones(self.n_compounds) / self.n_compounds
        
        # History for plotting
        self.P_history = []
        self.T_history = []
        self.F_history = []
        self.z_history = []

        # Ports this stream is connected to (for future extensibility)
        self.connected_ports = []
        self.arrayed_vars = {}  # Dict of {varname: [values]}
        self.printed_vars = {}  # Dict of {varname: [values]} (optional, for completeness)

    def update(self, P=None, T=None, F=None, z=None):
        if P is not None:
            self.P = P
        if T is not None:
            self.T = T
        if F is not None:
            self.F = F
        if z is not None:
            self.z = np.array(z) / np.sum(z)

    def store_history(self):
        self.P_history.append(self.P)
        self.T_history.append(self.T)
        self.F_history.append(self.F)
        # Store user-defined arrayed variables
        for var, value in self.arrayed_vars.items():
            hist_name = f"{var}_history"
            if not hasattr(self, hist_name):
                setattr(self, hist_name, [])
            getattr(self, hist_name).append(value)

    def run_user_code(self):
        # Try to load and execute user code for this stream
        code_file = 'user_code.py'
        code_block_name = f"# --- {self.name}_code ---"
        try:
            with open(code_file, 'r') as f:
                code_text = f.read()
            # Find the code block
            import re
            pattern = rf"{code_block_name}\n(.*?)(?=\n# ---|\Z)"
            match = re.search(pattern, code_text, re.DOTALL)
            if match:
                code = match.group(1).strip()
                local_vars = {
                    'self': self,
                    'CT': sim_settings.CT,
                    'DT': sim_settings.DT,
                    'RT': sim_settings.RT,
                    'TOTAL_TIME': sim_settings.TOTAL_TIME,
                }
                exec(code, {}, local_vars)
        except Exception as e:
            print(f"Error in stream user code for {self.name}: {e}")

    def connect_port(self, port):
        """Register a port (from a unit operation) that this stream is connected to."""
        self.connected_ports.append(port)

    def __repr__(self):
        return (f"<Stream {self.name}: P={self.P:.1f} Pa, T={self.T:.1f} K, "
                f"F={self.F:.3f} mol/s, z={self.z}>")
    @property
    def density(self):
        # Example for a mixture (ideal gas, modify for your needs)
        # self.z: mole fractions, self.compounds: list of names, self.T, self.P
        mix = Mixture(self.compounds, zs=self.z, T=self.T, P=self.P)
        return mix.rho  # kg/m3

    @property
    def enthalpy(self):
        mix = Mixture(self.compounds, zs=self.z, T=self.T, P=self.P)
        return mix.H  # J/mol or J/kg depending on thermo version

    @property
    def entropy(self):
        mix = Mixture(self.compounds, zs=self.z, T=self.T, P=self.P)
        return mix.S  # J/mol/K or J/kg/K