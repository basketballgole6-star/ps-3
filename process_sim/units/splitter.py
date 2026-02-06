# units/splitter.py

import numpy as np
from .unit_operation import UnitOperation

class Splitter(UnitOperation):
    """
    Splits the inlet flow into multiple outlet ports according to split fractions.
    """
    def __init__(self, name, n_outlets=2, split_fractions=None):
        super().__init__(name, n_inlets=1, n_outlets=n_outlets)
        # Default: equal split
        if split_fractions is None:
            split_fractions = [1.0 / n_outlets] * n_outlets
        split_fractions = np.array(split_fractions)
        split_fractions = split_fractions / np.sum(split_fractions)
        self.parameters = {'split_fractions': split_fractions}

    def calculate(self):
        inlet = self.inlet_ports[0].stream
        if inlet is None:
            return
        split_fractions = self.parameters['split_fractions']
        for i, port in enumerate(self.outlet_ports):
            outlet = port.stream
            if outlet is not None:
                F_out = inlet.F * split_fractions[i]
                outlet.update(F=F_out, P=inlet.P, T=inlet.T, z=inlet.z)
