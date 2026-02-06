import logging
import numpy as np

logger = logging.getLogger(__name__)

def get_constants_and_properties(components):
    """
    Get thermo constants and properties for a list of component names.
    Returns (constants, properties) tuple.
    """
    from thermo import ChemicalConstantsPackage
    return ChemicalConstantsPackage.from_IDs(components)

class Flasher:
    """
    Perform flash calculations for mixtures using thermo FlashVL.
    Supports multiple flash specifications: TP, PH, PS.
    """
    def __init__(self, components, composition, pressure=None, temperature=None, eos_choice='PRMIX'):
        from thermo import CEOSGas, CEOSLiquid, PRMIX, SRKMIX, RKMIX, APISRKMIX, TWUPRMIX, FlashVL
        from thermo.interaction_parameters import IPDB

        self.components = components
        self.composition = composition
        self.pressure = pressure
        self.temperature = temperature
        self.eos_choice = eos_choice

        # Validate composition
        total_composition = sum(composition.get(comp, 0) for comp in components)
        if not (0.99 <= total_composition <= 1.01):
            logger.warning(f"Composition sum is {total_composition}, should be ~1.0")

        # Convert composition to list
        self.zs = [composition.get(comp, 0) for comp in components]

        # Get constants and properties
        self.constants, self.properties = get_constants_and_properties(components)

        # Get interaction parameters (kij)
        kijs = None
        try:
            kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', self.constants.CASs, 'kij')
        except Exception:
            logger.debug("No kij parameters found, using None")

        # Setup EOS
        eos_map = {
            'PRMIX': PRMIX,
            'SRKMIX': SRKMIX,
            'RKMIX': RKMIX,
            'APISRKMIX': APISRKMIX,
            'TWUPRMIX': TWUPRMIX
        }
        eos_class = eos_map.get(eos_choice, PRMIX)
        eos_kwargs = {
            'Pcs': self.constants.Pcs,
            'Tcs': self.constants.Tcs,
            'omegas': self.constants.omegas,
            'kijs': kijs
        }

        # Create gas and liquid phases
        self.gas = CEOSGas(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases)
        self.liquid = CEOSLiquid(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases)

        # Create flasher
        self.flasher = FlashVL(self.constants, self.properties, liquid=self.liquid, gas=self.gas)

    def _format_result(self, flash_state, flash_type):
        try:
            vapor_fraction = getattr(flash_state, 'gas_beta', 1.0)
            liquid_fraction = getattr(flash_state, 'liquid_beta', 1.0 - vapor_fraction)
            if vapor_fraction >= 0.9999:
                phase = 'vapor'
            elif vapor_fraction <= 0.0001:
                phase = 'liquid'
            else:
                phase = 'two-phase'
            result = {
                'flash_type': flash_type,
                'vapor_fraction': vapor_fraction,
                'liquid_fraction': liquid_fraction,
                'temperature': flash_state.T,
                'pressure': flash_state.P,
                'phase': phase,
                'enthalpy': flash_state.H(),
                'entropy': flash_state.S(),
                'density': flash_state.rho() if hasattr(flash_state, 'rho') else None,
                'molar_volume': flash_state.V() if hasattr(flash_state, 'V') else None,
                'molar_weight': sum(self.constants.MWs[i] * self.zs[i] for i in range(len(self.zs))),
            }
            return result
        except Exception as e:
            logger.error(f"Error formatting result: {e}")
            return None

    def flash_tp(self, temperature, pressure):
        try:
            logger.debug(f"Performing TP-flash at T={temperature}K, P={pressure}Pa")
            flash_state = self.flasher.flash(T=temperature, P=pressure, zs=self.zs)
            return self._format_result(flash_state, 'TP')
        except Exception as e:
            logger.error(f"Error in TP-flash: {e}")
            return None

    def flash_ph(self, pressure, enthalpy):
        try:
            logger.debug(f"Performing PH-flash at P={pressure}Pa, H={enthalpy}J/mol")
            flash_state = self.flasher.flash(H=enthalpy, P=pressure, zs=self.zs)
            return self._format_result(flash_state, 'PH')
        except Exception as e:
            logger.error(f"Error in PH-flash: {e}")
            return None

    def flash_ps(self, pressure, entropy):
        try:
            logger.debug(f"Performing PS-flash at P={pressure}Pa, S={entropy}J/mol·K")
            flash_state = self.flasher.flash(S=entropy, P=pressure, zs=self.zs)
            return self._format_result(flash_state, 'PS')
        except Exception as e:
            logger.error(f"Error in PS-flash: {e}")
            return None

def create_flasher(components, composition, pressure=None, temperature=None, eos_choice='PRMIX'):
    return Flasher(components, composition, pressure, temperature, eos_choice)

def isenthalpic_expansion(components, composition, inlet_pressure, inlet_temperature, outlet_pressure, eos_choice='PRMIX'):
    try:
        logger.info(f"Isenthalpic expansion from {inlet_pressure/1e5:.2f} bar to {outlet_pressure/1e5:.2f} bar")
        zs = [composition.get(comp, 0) for comp in components]
        constants, properties = get_constants_and_properties(components)
        from thermo import CEOSGas, CEOSLiquid, PRMIX, SRKMIX, RKMIX, APISRKMIX, TWUPRMIX, FlashVL
        from thermo.interaction_parameters import IPDB
        eos_map = {
            'PRMIX': PRMIX,
            'SRKMIX': SRKMIX,
            'RKMIX': RKMIX,
            'APISRKMIX': APISRKMIX,
            'TWUPRMIX': TWUPRMIX
        }
        eos_class = eos_map.get(eos_choice, PRMIX)
        kijs = None
        try:
            kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        except Exception:
            pass
        eos_kwargs = {
            'Pcs': constants.Pcs,
            'Tcs': constants.Tcs,
            'omegas': constants.omegas,
            'kijs': kijs
        }
        gas = CEOSGas(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        inlet_state = flasher.flash(T=inlet_temperature, P=inlet_pressure, zs=zs)
        H_inlet = inlet_state.H()
        logger.debug(f"Inlet: T={inlet_temperature}K, P={inlet_pressure}Pa, H={H_inlet:.2f}J/mol")
        outlet_state = flasher.flash(H=H_inlet, P=outlet_pressure, zs=zs)
        vapor_fraction = getattr(outlet_state, 'gas_beta', 1.0)
        liquid_fraction = getattr(outlet_state, 'liquid_beta', 1.0 - vapor_fraction)
        if vapor_fraction >= 0.9999:
            phase = 'vapor'
        elif vapor_fraction <= 0.0001:
            phase = 'liquid'
        else:
            phase = 'two-phase'
        result = {
            'flash_type': 'PH',
            'temperature': outlet_state.T,
            'pressure': outlet_state.P,
            'enthalpy': outlet_state.H(),
            'entropy': outlet_state.S(),
            'vapor_fraction': vapor_fraction,
            'liquid_fraction': liquid_fraction,
            'phase': phase,
            'density': outlet_state.rho() if hasattr(outlet_state, 'rho') else None,
            'molar_volume': outlet_state.V() if hasattr(outlet_state, 'V') else None,
            'molar_weight': sum(constants.MWs[i] * zs[i] for i in range(len(zs))),
        }
        delta_T = inlet_temperature - result['temperature']
        logger.info(f"Outlet: T={result['temperature']:.2f}K, P={result['pressure']/1e5:.2f}bar, ΔT={delta_T:.2f}K")
        return result
    except Exception as e:
        logger.error(f"Error in isenthalpic expansion: {e}")
        import traceback
        traceback.print_exc()
        return None

def adiabatic_isentropic_expansion(components, composition, inlet_pressure, inlet_temperature, outlet_pressure, eos_choice='PRMIX'):
    try:
        logger.info(f"Isentropic expansion from {inlet_pressure/1e5:.2f} bar to {outlet_pressure/1e5:.2f} bar")
        zs = [composition.get(comp, 0) for comp in components]
        constants, properties = get_constants_and_properties(components)
        from thermo import CEOSGas, CEOSLiquid, PRMIX, SRKMIX, RKMIX, APISRKMIX, TWUPRMIX, FlashVL
        from thermo.interaction_parameters import IPDB
        eos_map = {
            'PRMIX': PRMIX,
            'SRKMIX': SRKMIX,
            'RKMIX': RKMIX,
            'APISRKMIX': APISRKMIX,
            'TWUPRMIX': TWUPRMIX
        }
        eos_class = eos_map.get(eos_choice, PRMIX)
        kijs = None
        try:
            kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        except Exception:
            pass
        eos_kwargs = {
            'Pcs': constants.Pcs,
            'Tcs': constants.Tcs,
            'omegas': constants.omegas,
            'kijs': kijs
        }
        gas = CEOSGas(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        inlet_state = flasher.flash(T=inlet_temperature, P=inlet_pressure, zs=zs)
        S_inlet = inlet_state.S()
        logger.debug(f"Inlet: T={inlet_temperature}K, P={inlet_pressure}Pa, S={S_inlet:.2f}J/mol·K")
        outlet_state = flasher.flash(S=S_inlet, P=outlet_pressure, zs=zs)
        vapor_fraction = getattr(outlet_state, 'gas_beta', 1.0)
        liquid_fraction = getattr(outlet_state, 'liquid_beta', 1.0 - vapor_fraction)
        if vapor_fraction >= 0.9999:
            phase = 'vapor'
        elif vapor_fraction <= 0.0001:
            phase = 'liquid'
        else:
            phase = 'two-phase'
        result = {
            'flash_type': 'PS',
            'temperature': outlet_state.T,
            'pressure': outlet_state.P,
            'enthalpy': outlet_state.H(),
            'entropy': outlet_state.S(),
            'vapor_fraction': vapor_fraction,
            'liquid_fraction': liquid_fraction,
            'phase': phase,
            'density': outlet_state.rho() if hasattr(outlet_state, 'rho') else None,
            'molar_volume': outlet_state.V() if hasattr(outlet_state, 'V') else None,
            'molar_weight': sum(constants.MWs[i] * zs[i] for i in range(len(zs))),
        }
        delta_T = inlet_temperature - result['temperature']
        logger.info(f"Outlet: T={result['temperature']:.2f}K, P={result['pressure']/1e5:.2f}bar, ΔT={delta_T:.2f}K")
        return result
    except Exception as e:
        logger.error(f"Error in isentropic expansion: {e}")
        import traceback
        traceback.print_exc()
        return None
