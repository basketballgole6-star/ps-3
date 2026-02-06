"""
Thermodynamic Property Extraction Module - FIXED VERSION
Uses thermo library for Peng-Robinson EOS based properties
All unit operations and streams use functions from this module

FIXES:
1. Robust error handling for missing properties
2. Better fallback calculations for Z-factor and Cp
3. Try-except blocks with more informative error messages
"""

from thermo import Chemical, Mixture
import numpy as np
import logging

# Setup logging
logger = logging.getLogger(__name__)


class PropertyExtractor:
    """
    Extract thermodynamic properties for mixtures using thermo library
    Peng-Robinson EOS based calculations
    """
    
    def __init__(self, components, composition, pressure, temperature):
        """
        Initialize property extractor for a mixture
        
        Args:
            components (list): List of component names (e.g., ['methane', 'ethane', 'propane'])
            composition (dict): Mole fractions {component: fraction}
            pressure (float): Pressure in Pa
            temperature (float): Temperature in K
        """
        self.components = components
        self.composition = composition
        self.pressure = pressure
        self.temperature = temperature
        
        # Validate composition
        total_composition = sum(composition.get(comp, 0) for comp in components)
        if not (0.99 <= total_composition <= 1.01):
            logger.warning(f"Composition sum is {total_composition}, should be ~1.0")
        
        try:
            # Create thermo Mixture object
            self.mixture = Mixture(
                IDs=components,
                zs=[composition.get(comp, 0) for comp in components],
                P=pressure,
                T=temperature
            )
            self.valid = True
        except Exception as e:
            logger.error(f"Error creating mixture: {e}")
            self.valid = False
            self.mixture = None
    
    # === BASIC PROPERTIES ===
    
    def get_density(self):
        """Get mixture density (kg/m³)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'rho') and self.mixture.rho is not None:
                    return self.mixture.rho
                # Fallback calculation
                elif hasattr(self.mixture, 'MW') and hasattr(self.mixture, 'Vm'):
                    return self.mixture.MW / self.mixture.Vm if self.mixture.Vm > 0 else None
            return None
        except Exception as e:
            logger.warning(f"Error getting density: {e}")
            return None
    
    def get_molar_volume(self):
        """Get mixture molar volume (m³/kmol)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'Vm'):
                return self.mixture.Vm
            return None
        except Exception as e:
            logger.warning(f"Error getting molar volume: {e}")
            return None
    
    def get_molar_weight(self):
        """Get mixture molar weight (kg/kmol)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'MW'):
                return self.mixture.MW
            # Fallback: calculate from components
            MW = sum(self.composition.get(comp, 0) * Chemical(comp).MW 
                    for comp in self.components)
            return MW if MW > 0 else None
        except Exception as e:
            logger.warning(f"Error getting molar weight: {e}")
            return None
    
    def get_critical_temp(self):
        """Get mixture critical temperature (K) using Kay's rule"""
        try:
            Tc = sum(self.composition.get(comp, 0) * Chemical(comp).Tc 
                    for comp in self.components)
            return Tc
        except Exception as e:
            logger.warning(f"Error getting critical temperature: {e}")
            return None
    
    def get_critical_pressure(self):
        """Get mixture critical pressure (Pa) using Kay's rule"""
        try:
            Pc = sum(self.composition.get(comp, 0) * Chemical(comp).Pc 
                    for comp in self.components)
            return Pc
        except Exception as e:
            logger.warning(f"Error getting critical pressure: {e}")
            return None
    
    def get_critical_volume(self):
        """Get mixture critical volume (m³/kmol)"""
        try:
            Vc = sum(self.composition.get(comp, 0) * Chemical(comp).Vc 
                    for comp in self.components)
            return Vc
        except Exception as e:
            logger.warning(f"Error getting critical volume: {e}")
            return None
    
    def get_acentric_factor(self):
        """Get mixture acentric factor (dimensionless)"""
        try:
            omega = sum(self.composition.get(comp, 0) * Chemical(comp).omega 
                       for comp in self.components)
            return omega
        except Exception as e:
            logger.warning(f"Error getting acentric factor: {e}")
            return None
    
    # === ENTHALPY & ENTROPY ===
    
    def get_enthalpy(self):
        """Get mixture enthalpy (J/kmol)"""
        try:
            if self.valid and self.mixture is not None:
                # Try different methods to get enthalpy
                if hasattr(self.mixture, 'H') and self.mixture.H is not None:
                    return self.mixture.H
                elif hasattr(self.mixture, 'Hm') and self.mixture.Hm is not None:
                    return self.mixture.Hm * 1000  # Convert J/mol to J/kmol
            return None
        except Exception as e:
            logger.warning(f"Error getting enthalpy: {e}")
            return None
    
    def get_enthalpy_mass(self):
        """Get mixture specific enthalpy (J/kg)"""
        try:
            H = self.get_enthalpy()
            MW = self.get_molar_weight()
            if H is not None and MW is not None and MW > 0:
                return H / MW
            return None
        except Exception as e:
            logger.warning(f"Error getting mass enthalpy: {e}")
            return None
    
    def get_entropy(self):
        """Get mixture entropy (J/kmol·K)"""
        try:
            if self.valid and self.mixture is not None:
                # Try different methods
                if hasattr(self.mixture, 'S') and self.mixture.S is not None:
                    return self.mixture.S
                elif hasattr(self.mixture, 'Sm') and self.mixture.Sm is not None:
                    return self.mixture.Sm * 1000  # Convert J/mol·K to J/kmol·K
            return None
        except Exception as e:
            logger.warning(f"Error getting entropy: {e}")
            return None
    
    def get_entropy_mass(self):
        """Get mixture specific entropy (J/kg·K)"""
        try:
            S = self.get_entropy()
            MW = self.get_molar_weight()
            if S is not None and MW is not None and MW > 0:
                return S / MW
            return None
        except Exception as e:
            logger.warning(f"Error getting mass entropy: {e}")
            return None
    
    def get_gibbs_energy(self):
        """Get mixture Gibbs energy (J/kmol)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'G') and self.mixture.G is not None:
                    return self.mixture.G
                elif hasattr(self.mixture, 'Gm') and self.mixture.Gm is not None:
                    return self.mixture.Gm * 1000
            return None
        except Exception as e:
            logger.warning(f"Error getting Gibbs energy: {e}")
            return None
    
    def get_internal_energy(self):
        """Get mixture internal energy (J/kmol)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'U') and self.mixture.U is not None:
                    return self.mixture.U
                elif hasattr(self.mixture, 'Um') and self.mixture.Um is not None:
                    return self.mixture.Um * 1000
            return None
        except Exception as e:
            logger.warning(f"Error getting internal energy: {e}")
            return None
    
    # === HEAT CAPACITY ===
    
    def get_cp(self):
        """Get mixture heat capacity at constant pressure (J/kmol·K)"""
        try:
            if self.valid and self.mixture is not None:
                # Try multiple approaches
                if hasattr(self.mixture, 'Cp') and self.mixture.Cp is not None:
                    return self.mixture.Cp
                elif hasattr(self.mixture, 'Cpm') and self.mixture.Cpm is not None:
                    return self.mixture.Cpm * 1000  # J/mol·K to J/kmol·K
                elif hasattr(self.mixture, 'Cpg') and self.mixture.Cpg is not None:
                    return self.mixture.Cpg
                elif hasattr(self.mixture, 'Cpl') and self.mixture.Cpl is not None:
                    return self.mixture.Cpl
                
                # Fallback: Calculate from components
                logger.debug("Using component-based Cp calculation")
                Cp_mix = 0.0
                for i, comp in enumerate(self.components):
                    try:
                        chem = Chemical(comp, P=self.pressure, T=self.temperature)
                        comp_cp = None
                        if hasattr(chem, 'Cp') and chem.Cp is not None:
                            comp_cp = chem.Cp
                        elif hasattr(chem, 'Cpm') and chem.Cpm is not None:
                            comp_cp = chem.Cpm * 1000
                        elif hasattr(chem, 'Cpg') and chem.Cpg is not None:
                            comp_cp = chem.Cpg
                        
                        if comp_cp is not None:
                            Cp_mix += self.composition.get(comp, 0) * comp_cp
                    except Exception as e:
                        logger.debug(f"Could not get Cp for {comp}: {e}")
                        continue
                
                if Cp_mix > 0:
                    return Cp_mix
            
            # Ultimate fallback: estimate from ideal gas Cp
            logger.debug("Using ideal gas Cp estimate")
            R = 8.314  # J/mol·K
            return 3.5 * R * 1000  # Rough estimate for diatomic/polyatomic gas (J/kmol·K)
            
        except Exception as e:
            logger.warning(f"Error getting Cp: {e}")
            # Return a reasonable default
            return 35000.0  # ~3.5R in J/kmol·K
    
    def get_cp_mass(self):
        """Get mixture specific heat capacity at constant pressure (J/kg·K)"""
        try:
            Cp = self.get_cp()
            MW = self.get_molar_weight()
            if Cp is not None and MW is not None and MW > 0:
                return Cp / MW
            return None
        except Exception as e:
            logger.warning(f"Error getting mass Cp: {e}")
            return None
    
    def get_cv(self):
        """Get mixture heat capacity at constant volume (J/kmol·K)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'Cv') and self.mixture.Cv is not None:
                    return self.mixture.Cv
                elif hasattr(self.mixture, 'Cvm') and self.mixture.Cvm is not None:
                    return self.mixture.Cvm * 1000
            
            # Fallback: Cv = Cp - R (for ideal gas)
            Cp = self.get_cp()
            if Cp is not None:
                R = 8.314 * 1000  # J/kmol·K
                return Cp - R
            
            return None
        except Exception as e:
            logger.warning(f"Error getting Cv: {e}")
            return None
    
    def get_cv_mass(self):
        """Get mixture specific heat capacity at constant volume (J/kg·K)"""
        try:
            Cv = self.get_cv()
            MW = self.get_molar_weight()
            if Cv is not None and MW is not None and MW > 0:
                return Cv / MW
            return None
        except Exception as e:
            logger.warning(f"Error getting mass Cv: {e}")
            return None
    
    def get_gamma(self):
        """Get mixture heat capacity ratio (Cp/Cv) - dimensionless"""
        try:
            cp = self.get_cp()
            cv = self.get_cv()
            if cp and cv and cv > 0:
                return cp / cv
            # Fallback: typical value for gases
            return 1.4
        except Exception as e:
            logger.warning(f"Error getting gamma: {e}")
            return 1.4  # Reasonable default for diatomic gas
    
    # === VISCOSITY & THERMAL PROPERTIES ===
    
    def get_viscosity(self):
        """Get mixture viscosity (Pa·s)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'mu') and self.mixture.mu is not None:
                    return self.mixture.mu
                elif hasattr(self.mixture, 'mug') and self.mixture.mug is not None:
                    return self.mixture.mug
                elif hasattr(self.mixture, 'mul') and self.mixture.mul is not None:
                    return self.mixture.mul
            return None
        except Exception as e:
            logger.warning(f"Error getting viscosity: {e}")
            return None
    
    def get_thermal_conductivity(self):
        """Get mixture thermal conductivity (W/m·K)"""
        try:
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'k') and self.mixture.k is not None:
                    return self.mixture.k
                elif hasattr(self.mixture, 'kg') and self.mixture.kg is not None:
                    return self.mixture.kg
                elif hasattr(self.mixture, 'kl') and self.mixture.kl is not None:
                    return self.mixture.kl
            return None
        except Exception as e:
            logger.warning(f"Error getting thermal conductivity: {e}")
            return None
    
    def get_surface_tension(self):
        """Get surface tension (N/m)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'sigma'):
                return self.mixture.sigma
            return None
        except Exception as e:
            logger.warning(f"Error getting surface tension: {e}")
            return None
    
    # === FUGACITY & COEFFICIENTS ===
    
    def get_fugacity(self):
        """Get component fugacities in mixture (Pa)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'fugacities'):
                return self.mixture.fugacities
            return None
        except Exception as e:
            logger.warning(f"Error getting fugacities: {e}")
            return None
    
    def get_fugacity_coefficient(self):
        """Get component fugacity coefficients (dimensionless)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'phi'):
                return self.mixture.phi
            return None
        except Exception as e:
            logger.warning(f"Error getting fugacity coefficients: {e}")
            return None
    
    def get_activity_coefficient(self):
        """Get component activity coefficients (for liquid phase)"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'gammas'):
                return self.mixture.gammas
            return None
        except Exception as e:
            logger.warning(f"Error getting activity coefficients: {e}")
            return None
    
    # === COMPRESSIBILITY & FLOW PROPERTIES ===
    
    def get_compressibility_factor(self):
        """
        Get compressibility factor (Z-factor) - IMPROVED VERSION
        Measures deviation from ideal gas behavior
        Z = PV / nRT
        Z = 1 for ideal gas
        Z < 1 for attractive forces dominant
        Z > 1 for repulsive forces dominant
        
        Returns: float (defaults to 1.0 if unavailable)
        """
        try:
            # Method 1: Direct Z property
            if self.valid and self.mixture is not None:
                if hasattr(self.mixture, 'Z') and self.mixture.Z is not None:
                    return self.mixture.Z
                
                # Method 2: Calculate from PV = ZnRT
                if hasattr(self.mixture, 'Vm') and self.mixture.Vm is not None:
                    R = 8.314  # J/mol·K
                    # Z = PV/(RT) where V is molar volume
                    Z = (self.pressure * self.mixture.Vm) / (R * self.temperature)
                    if 0.1 <= Z <= 10:  # Sanity check
                        return Z
                
                # Method 3: Use virial coefficient if available
                if hasattr(self.mixture, 'Bm') and self.mixture.Bm is not None:
                    R = 8.314  # J/mol·K
                    B = self.mixture.Bm  # m³/mol
                    Z = 1 + (self.pressure * B) / (R * self.temperature)
                    if 0.1 <= Z <= 10:  # Sanity check
                        return Z
            
            # Method 4: Estimate from reduced properties (Pitzer correlation)
            Tc = self.get_critical_temp()
            Pc = self.get_critical_pressure()
            omega = self.get_acentric_factor()
            
            if Tc and Pc and self.temperature > 0 and Pc > 0:
                Tr = self.temperature / Tc  # Reduced temperature
                Pr = self.pressure / Pc      # Reduced pressure
                
                # Simple Pitzer correlation for Z
                if omega is not None:
                    Z0 = 1 + (0.083 - 0.422 / (Tr ** 1.6)) * Pr
                    Z1 = (0.139 - 0.172 / (Tr ** 4.2)) * Pr
                    Z = Z0 + omega * Z1
                    if 0.1 <= Z <= 10:  # Sanity check
                        return Z
            
            # Default: ideal gas
            logger.debug("Using ideal gas assumption (Z=1.0)")
            return 1.0
            
        except Exception as e:
            logger.warning(f"Error getting compressibility factor: {e}")
            return 1.0  # Safe default
    
    def get_speed_of_sound(self):
        """
        Get speed of sound in mixture (m/s)
        Important for choked flow calculations
        c = sqrt(γ * R * T)
        
        Returns: float or None
        """
        try:
            if self.valid and self.mixture is not None:
                # Method 1: Direct property
                if hasattr(self.mixture, 'speed_of_sound') and self.mixture.speed_of_sound is not None:
                    return self.mixture.speed_of_sound
            
            # Method 2: Calculate from gamma
            gamma = self.get_gamma()
            MW = self.get_molar_weight()
            
            if gamma is not None and MW is not None and MW > 0:
                R = 8.314 / MW * 1000  # J/kg·K (specific gas constant)
                c = np.sqrt(gamma * R * self.temperature)
                return c
            
            return None
        except Exception as e:
            logger.warning(f"Error getting speed of sound: {e}")
            return None
    
    def get_reduced_pressure(self):
        """Get reduced pressure (Pr = P / Pc)"""
        try:
            Pc = self.get_critical_pressure()
            if Pc and Pc > 0:
                return self.pressure / Pc
            return None
        except Exception as e:
            logger.warning(f"Error getting reduced pressure: {e}")
            return None
    
    def get_reduced_temperature(self):
        """Get reduced temperature (Tr = T / Tc)"""
        try:
            Tc = self.get_critical_temp()
            if Tc and Tc > 0:
                return self.temperature / Tc
            return None
        except Exception as e:
            logger.warning(f"Error getting reduced temperature: {e}")
            return None
    
    # === PHASE IDENTIFICATION ===
    
    def get_phase_state(self):
        """Determine phase state of mixture"""
        try:
            if self.valid and self.mixture is not None and hasattr(self.mixture, 'phase'):
                phase = self.mixture.phase
                if phase in ['g', 'vapor', 'V', 'gas']:
                    return 'vapor'
                elif phase in ['l', 'liquid', 'L']:
                    return 'liquid'
                elif phase in ['l/g', 'two-phase', 'two phase']:
                    return 'two-phase'
                return phase
            return 'unknown'
        except Exception as e:
            logger.warning(f"Error getting phase state: {e}")
            return 'unknown'
    
    def is_two_phase(self):
        """Check if mixture is in two-phase region"""
        phase = self.get_phase_state()
        return phase == 'two-phase'
    
    def is_vapor(self):
        """Check if mixture is vapor phase"""
        phase = self.get_phase_state()
        return phase == 'vapor'
    
    def is_liquid(self):
        """Check if mixture is liquid phase"""
        phase = self.get_phase_state()
        return phase == 'liquid'
    
    # === COMPONENT PROPERTIES ===
    
    def get_component_properties(self, component):
        """Get individual component properties"""
        try:
            chem = Chemical(component, P=self.pressure, T=self.temperature)
            return {
                'name': component,
                'MW': chem.MW if hasattr(chem, 'MW') else None,
                'density': chem.rho if hasattr(chem, 'rho') else None,
                'enthalpy': chem.H if hasattr(chem, 'H') else None,
                'entropy': chem.S if hasattr(chem, 'S') else None,
                'Cp': chem.Cp if hasattr(chem, 'Cp') else None,
                'Cv': chem.Cv if hasattr(chem, 'Cv') else None,
                'viscosity': chem.mu if hasattr(chem, 'mu') else None,
                'thermal_conductivity': chem.k if hasattr(chem, 'k') else None,
                'Tc': chem.Tc if hasattr(chem, 'Tc') else None,
                'Pc': chem.Pc if hasattr(chem, 'Pc') else None,
                'Vc': chem.Vc if hasattr(chem, 'Vc') else None,
                'omega': chem.omega if hasattr(chem, 'omega') else None,
            }
        except Exception as e:
            logger.warning(f"Error getting component {component} properties: {e}")
            return None
    
    def get_all_components_properties(self):
        """Get properties for all components in the mixture"""
        all_props = {}
        for component in self.components:
            all_props[component] = self.get_component_properties(component)
        return all_props


# === CONVENIENCE FUNCTIONS ===

def create_property_extractor(components, composition, pressure, temperature):
    """Create and return a PropertyExtractor object"""
    return PropertyExtractor(components, composition, pressure, temperature)


def get_all_properties(components, composition, pressure, temperature):
    """Get all mixture properties at once"""
    props = PropertyExtractor(components, composition, pressure, temperature)
    
    return {
        'density': props.get_density(),
        'molar_volume': props.get_molar_volume(),
        'molar_weight': props.get_molar_weight(),
        'temperature': temperature,
        'pressure': pressure,
        'phase': props.get_phase_state(),
        'enthalpy': props.get_enthalpy(),
        'entropy': props.get_entropy(),
        'gibbs_energy': props.get_gibbs_energy(),
        'internal_energy': props.get_internal_energy(),
        'Cp': props.get_cp(),
        'Cv': props.get_cv(),
        'gamma': props.get_gamma(),
        'viscosity': props.get_viscosity(),
        'thermal_conductivity': props.get_thermal_conductivity(),
        'surface_tension': props.get_surface_tension(),
        'fugacity': props.get_fugacity(),
        'fugacity_coefficient': props.get_fugacity_coefficient(),
        'activity_coefficient': props.get_activity_coefficient(),
        'critical_temperature': props.get_critical_temp(),
        'critical_pressure': props.get_critical_pressure(),
        'critical_volume': props.get_critical_volume(),
        'acentric_factor': props.get_acentric_factor(),
        'compressibility_factor': props.get_compressibility_factor(),
        'speed_of_sound': props.get_speed_of_sound(),
        'reduced_pressure': props.get_reduced_pressure(),
        'reduced_temperature': props.get_reduced_temperature(),
    }