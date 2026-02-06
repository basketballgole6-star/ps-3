"""
SIMPLE FLASH TEST - One case only
"""

from thermoprops.flasher import isenthalpic_expansion

print("="*60)
print("SIMPLE ISENTHALPIC EXPANSION TEST")
print("="*60)

# Test conditions
components = ['methane', 'ethane', 'propane']
composition = {'methane': 0.5, 'ethane': 0.3, 'propane': 0.2}

P_in = 10e5   # 10 bar
T_in = 300    # 300 K
P_out = 7e5   # 7 bar

print(f"\nInput:")
print(f"  Pressure:    {P_in/1e5:.1f} bar")
print(f"  Temperature: {T_in:.1f} K")
print(f"  Composition: CH4={composition['methane']}, C2H6={composition['ethane']}, C3H8={composition['propane']}")

print(f"\nExpansion to: {P_out/1e5:.1f} bar")

# Run flash
result = isenthalpic_expansion(
    components=components,
    composition=composition,
    inlet_pressure=P_in,
    inlet_temperature=T_in,
    outlet_pressure=P_out
)

print("\n" + "="*60)
print("RESULT:")
print("="*60)

if result is None:
    print("❌ FLASH FAILED - returned None")
else:
    T_out = result.get('temperature')
    P_out_result = result.get('pressure')
    delta_T = T_in - T_out
    
    print(f"  Temperature:     {T_out:.2f} K")
    print(f"  Pressure:        {P_out_result/1e5:.2f} bar")
    print(f"  ΔT:              {delta_T:.2f} K")
    print(f"  Vapor fraction:  {result.get('vapor_fraction', 'N/A'):.4f}")
    print(f"  Phase:           {result.get('phase', 'N/A')}")
    
    if abs(delta_T) > 0.1:
        print(f"\n✓ SUCCESS: Temperature dropped by {delta_T:.2f} K!")
    else:
        print(f"\n⚠️  Temperature change is minimal ({delta_T:.2f} K)")
        print("   This may be normal for these conditions (near-ideal gas)")

print("="*60)