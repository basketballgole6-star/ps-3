from thermoprops.flasher import Flasher,isenthalpic_expansion,adiabatic_isentropic_expansion
from thermoprops.property_extractor import PropertyExtractor, get_all_properties

# Example 1: Get mixture properties
components = ['methane', 'ethane', 'propane']
composition = {'methane': 0.7, 'ethane': 0.2, 'propane': 0.1}
pressure = 10e5  # 10 bar in Pa
temperature = 300  # K

props = PropertyExtractor(components, composition, pressure, temperature)
print(f"Density: {props.get_density()} kg/m³")
print(f"Enthalpy: {props.get_enthalpy()} J/kmol")
print(f"Entropy: {props.get_entropy()} J/kmol·K")

# Example 2: Get all properties at once
all_props = get_all_properties(components, composition, pressure, temperature)
print(all_props)

# Example 3: TP-flash
flasher = Flasher(components, composition)
flash_result = flasher.flash_tp(temperature=300, pressure=5e5)
print(f"Vapor fraction: {flash_result['vapor_fraction']}")
print(f"Temperature: {flash_result['temperature']} K")

# Example 4: Isenthalpic expansion (valve, throttling)
valve_outlet = isenthalpic_expansion(
    components, 
    composition,
    inlet_pressure=10e5,
    inlet_temperature=300,
    outlet_pressure=5e5
)
print(f"After valve - Vapor fraction: {valve_outlet['vapor_fraction']}")
print(f"After valve - Temperature: {valve_outlet['temperature']} K")

# Example 5: Isentropic expansion (turbine)
turbine_outlet = adiabatic_isentropic_expansion(
    components,
    composition,
    inlet_pressure=10e5,
    inlet_temperature=300,
    outlet_pressure=5e5
)
print(f"After turbine - Temperature: {turbine_outlet['temperature']} K")
print(f"After turbine - Enthalpy: {turbine_outlet['enthalpy']} J/kmol")