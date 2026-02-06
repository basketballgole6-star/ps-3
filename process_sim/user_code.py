"""
User Code for Process Simulation
Executed at specific points in unit operations
Define code blocks like:
# --- UnitName_before ---
# --- UnitName_after ---
"""

# ============================================================================
# FEED BOUNDARY - Source
# ============================================================================
# --- Feed_before ---
# Nothing to do before feed
# --- Feed_before_end ---

# --- Feed_after ---
self.arrayed_vars['P_out'] = self.outlet_ports[0].stream.P / 1e5
self.arrayed_vars['T_out'] = self.outlet_ports[0].stream.T
self.arrayed_vars['F_out'] = self.outlet_ports[0].stream.F
# --- Feed_after_end ---


# ============================================================================
# VALVE1 - Expansion Valve (10 bar → 8 bar)
# ============================================================================
# --- Valve1_before ---
inlet = self.inlet_ports[0].stream
self.arrayed_vars['P_in'] = inlet.P / 1e5
self.arrayed_vars['T_in'] = inlet.T
self.arrayed_vars['F_in'] = inlet.F
if CT>0:
    self.parameters['opening'] = max(1 - CT/10 ,0.0001)
else:
    self.parameters['opening'] = 1
self.arrayed_vars['opening'] = self.parameters['opening']
# --- Valve1_before_end ---

# --- Valve1_after ---
inlet = self.inlet_ports[0].stream
outlet = self.outlet_ports[0].stream

# Outlet conditions
self.arrayed_vars['P_out'] = outlet.P / 1e5
self.arrayed_vars['T_out'] = outlet.T
self.arrayed_vars['F_out'] = outlet.F
self.arrayed_vars['dP'] = self.results.get('dP', 0) / 1e5

# Thermodynamic properties
self.arrayed_vars['vapor_fraction'] = self.results.get('vapor_fraction', 0)
self.arrayed_vars['gamma'] = self.results.get('gamma', 1.4)
self.arrayed_vars['Z_inlet'] = self.results.get('Z_inlet', 1.0)
self.arrayed_vars['Z_outlet'] = self.results.get('Z_outlet', 1.0)
self.arrayed_vars['Cv_effective'] = self.results.get('Cv_effective', 0)

# Flash results
flash_result = self.results.get('flash_result', {})
self.arrayed_vars['H_outlet'] = flash_result.get('enthalpy', 0)
self.arrayed_vars['S_outlet'] = flash_result.get('entropy', 0)
self.arrayed_vars['rho_outlet'] = flash_result.get('density', 0)

# Flow characteristics
self.arrayed_vars['is_choked'] = int(self.results.get('is_choked', False))

# Status
self.arrayed_vars['status'] = 1 if self.results.get('status') == 'open' else 0

# Print summary
print(f"\n{'='*80}")
print(f"VALVE: {self.name}")
print(f"{'='*80}")
print(f"TIME: {CT:.4f}s")
print(f"\n{'INLET':-^80}")
print(f"  Pressure:    {inlet.P/1e5:>10.3f} bar")
print(f"  Temperature: {inlet.T:>10.2f} K")
print(f"  Flow Rate:   {inlet.F:>10.4f} kg/s")
print(f"\n{'VALVE SETTINGS':-^80}")
print(f"  Opening:     {self.parameters['opening']*100:>10.1f}%")
print(f"  Cv (base):   {self.parameters['Cv']:>10.2f}")
print(f"  Cv (eff):    {self.results.get('Cv_effective', 0):>10.4f}")
print(f"\n{'OUTLET':-^80}")
print(f"  Pressure:    {outlet.P/1e5:>10.3f} bar (target: 8.00)")
print(f"  Temperature: {outlet.T:>10.2f} K")
print(f"  Flow Rate:   {outlet.F:>10.4f} kg/s")
print(f"  Vapor Frac:  {self.results.get('vapor_fraction', 0):>10.2%}")
print(f"\n{'THERMODYNAMICS':-^80}")
print(f"  ΔP:          {self.results.get('dP', 0)/1e5:>10.3f} bar")
print(f"  ΔT:          {inlet.T - outlet.T:>10.2f} K")
print(f"  Z (inlet):   {self.results.get('Z_inlet', 1.0):>10.4f}")
print(f"  Z (outlet):  {self.results.get('Z_outlet', 1.0):>10.4f}")
print(f"  γ:           {self.results.get('gamma', 1.4):>10.4f}")
print(f"  Choked:      {'YES' if self.results.get('is_choked') else 'NO':>10}")
print(f"{'='*80}\n")
# --- Valve1_after_end ---


# ============================================================================
# VOLUME - Pressure Vessel/Accumulator
# ============================================================================
# --- Volume_1_before ---
inlet = self.inlet_ports[0].stream
self.arrayed_vars['inlet_P'] = inlet.P / 1e5
self.arrayed_vars['inlet_T'] = inlet.T
self.arrayed_vars['inlet_F'] = inlet.F
# --- Volume_1_before_end ---

# --- Volume_1_after ---
inlet = self.inlet_ports[0].stream
outlet = self.outlet_ports[0].stream

self.arrayed_vars['P'] = self.P / 1e5
self.arrayed_vars['T'] = self.T
self.arrayed_vars['F_in'] = inlet.F
self.arrayed_vars['F_out'] = outlet.F
self.arrayed_vars['dP_inlet_to_vessel'] = (inlet.P - self.P) / 1e5

print(f"\n{'='*80}")
print(f"VOLUME: {self.name}")
print(f"{'='*80}")
print(f"TIME: {CT:.4f}s")
print(f"\n{'INLET CONDITIONS':-^80}")
print(f"  Pressure:    {inlet.P/1e5:>10.3f} bar")
print(f"  Temperature: {inlet.T:>10.2f} K")
print(f"  Flow Rate:   {inlet.F:>10.4f} kg/s")
print(f"\n{'VESSEL STATE':-^80}")
print(f"  Volume:      {self.V_total:>10.2f} m³")
print(f"  Pressure:    {self.P/1e5:>10.3f} bar")
print(f"  Temperature: {self.T:>10.2f} K")
print(f"  Inlet-Vessel ΔP: {(inlet.P - self.P)/1e5:>10.3f} bar")
print(f"\n{'OUTLET CONDITIONS':-^80}")
print(f"  Pressure:    {outlet.P/1e5:>10.3f} bar")
print(f"  Temperature: {outlet.T:>10.2f} K")
print(f"  Flow Rate:   {outlet.F:>10.4f} kg/s")
print(f"{'='*80}\n")
# --- Volume_1_after_end ---


# ============================================================================
# VALVE2 - Relief Valve (8 bar → 1.5 bar)
# ============================================================================
# --- Valve2_before ---
inlet = self.inlet_ports[0].stream
self.arrayed_vars['P_in'] = inlet.P / 1e5
self.arrayed_vars['T_in'] = inlet.T
self.arrayed_vars['F_in'] = inlet.F
self.arrayed_vars['opening'] = self.parameters['opening']
# --- Valve2_before_end ---

# --- Valve2_after ---
inlet = self.inlet_ports[0].stream
outlet = self.outlet_ports[0].stream

# Outlet conditions
self.arrayed_vars['P_out'] = outlet.P / 1e5
self.arrayed_vars['T_out'] = outlet.T
self.arrayed_vars['F_out'] = outlet.F
self.arrayed_vars['dP'] = self.results.get('dP', 0) / 1e5

# Thermodynamic properties
self.arrayed_vars['vapor_fraction'] = self.results.get('vapor_fraction', 0)
self.arrayed_vars['gamma'] = self.results.get('gamma', 1.4)
self.arrayed_vars['Z_inlet'] = self.results.get('Z_inlet', 1.0)
self.arrayed_vars['Z_outlet'] = self.results.get('Z_outlet', 1.0)
self.arrayed_vars['Cv_effective'] = self.results.get('Cv_effective', 0)

# Flash results
flash_result = self.results.get('flash_result', {})
self.arrayed_vars['H_outlet'] = flash_result.get('enthalpy', 0)
self.arrayed_vars['S_outlet'] = flash_result.get('entropy', 0)
self.arrayed_vars['rho_outlet'] = flash_result.get('density', 0)

# Flow characteristics
self.arrayed_vars['is_choked'] = int(self.results.get('is_choked', False))

# Status
self.arrayed_vars['status'] = 1 if self.results.get('status') == 'open' else 0

# Print summary
print(f"\n{'='*80}")
print(f"VALVE: {self.name}")
print(f"{'='*80}")
print(f"TIME: {CT:.4f}s")
print(f"\n{'INLET':-^80}")
print(f"  Pressure:    {inlet.P/1e5:>10.3f} bar")
print(f"  Temperature: {inlet.T:>10.2f} K")
print(f"  Flow Rate:   {inlet.F:>10.4f} kg/s")
print(f"\n{'VALVE SETTINGS':-^80}")
print(f"  Opening:     {self.parameters['opening']*100:>10.1f}%")
print(f"  Cv (base):   {self.parameters['Cv']:>10.2f}")
print(f"  Cv (eff):    {self.results.get('Cv_effective', 0):>10.4f}")
print(f"\n{'OUTLET':-^80}")
print(f"  Pressure:    {outlet.P/1e5:>10.3f} bar (target: 1.50)")
print(f"  Temperature: {outlet.T:>10.2f} K")
print(f"  Flow Rate:   {outlet.F:>10.4f} kg/s")
print(f"  Vapor Frac:  {self.results.get('vapor_fraction', 0):>10.2%}")
print(f"\n{'THERMODYNAMICS':-^80}")
print(f"  ΔP:          {self.results.get('dP', 0)/1e5:>10.3f} bar")
print(f"  ΔT:          {inlet.T - outlet.T:>10.2f} K")
print(f"  Z (inlet):   {self.results.get('Z_inlet', 1.0):>10.4f}")
print(f"  Z (outlet):  {self.results.get('Z_outlet', 1.0):>10.4f}")
print(f"  γ:           {self.results.get('gamma', 1.4):>10.4f}")
print(f"  Choked:      {'YES' if self.results.get('is_choked') else 'NO':>10}")
print(f"{'='*80}\n")
# --- Valve2_after_end ---


# ============================================================================
# PRODUCT BOUNDARY - Exit
# ============================================================================
# --- Product_before ---
# Nothing to do before product
# --- Product_before_end ---

# --- Product_after ---
inlet = self.inlet_ports[0].stream
self.arrayed_vars['P_in'] = inlet.P / 1e5
self.arrayed_vars['T_in'] = inlet.T
self.arrayed_vars['F_in'] = inlet.F
# --- Product_after_end ---