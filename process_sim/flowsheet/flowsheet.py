"""
Improved Test Flowsheet: Natural Gas Expansion System
Format: Using FeedBoundary, Boundary, Valve, and Volume units

IMPROVEMENTS:
1. More moderate pressure drops (better for numerical stability)
2. Realistic valve Cv values
3. Better initial conditions
4. Added intermediate pressure levels

System Layout (IMPROVED):
┌──────────┐     ┌────────┐     ┌─────────┐     ┌────────┐     ┌──────────┐
│   Feed   │────▶│ Valve1 │────▶│ Volume  │────▶│ Valve2 │────▶│ Product  │
│10b,300K  │     │10→7bar │     │  7bar   │     │7→3bar  │     │3b,~295K  │
└──────────┘     └────────┘     └─────────┘     └────────┘     └──────────┘

KEY CHANGES FROM ORIGINAL:
- Valve1: 10→7 bar (was 10→8, more gradual drop)
- Volume: 7 bar setpoint (was 8 bar)
- Valve2: 7→3 bar (was 8→1.5, less extreme)
- Product: 3 bar boundary (was 1.5 bar)
- Larger Cv values for smoother flow
- Expected temp drop: 2-3K per valve (total ~5-6K)
"""

from units.stream import Stream
from units.feed import FeedBoundary
from units.valve import Valve
from units.volume import Volume
from units.boundary import Boundary


def build_flowsheet(sim, compounds):
    """
    Build the improved natural gas expansion test system
    """
    
    # =========================================================================
    # CREATE STREAMS WITH BETTER INITIAL CONDITIONS
    # =========================================================================
    
    # Stream 1: Feed outlet (10 bar, 300K)
    s1 = Stream(
        name='s1_FeedOutlet',
        compounds=compounds,
        P_init=10e5,           # 10 bar
        T_init=300,            # 300 K
        z_init=[0.5, 0.3, 0.2],  # CH4: 50%, C2H6: 30%, C3H8: 20%
        F_init=10.0            # 10 kg/s
    )
    
    # Stream 2: Valve1 outlet (7 bar, ~298K expected)
    s2 = Stream(
        name='s2_Valve1Outlet',
        compounds=compounds,
        P_init=7e5,            # 7 bar (more moderate drop)
        T_init=298,            # Slight temp drop expected
        z_init=[0.5, 0.3, 0.2],
        F_init=10.0            # Initialize with flow
    )
    
    # Stream 3: Volume outlet (7 bar, stabilized)
    s3 = Stream(
        name='s3_VolumeOutlet',
        compounds=compounds,
        P_init=7e5,            # 7 bar
        T_init=298,
        z_init=[0.5, 0.3, 0.2],
        F_init=10.0
    )
    
    # Stream 4: Valve2 outlet (3 bar, ~295K expected)
    s4 = Stream(
        name='s4_Valve2Outlet',
        compounds=compounds,
        P_init=3e5,            # 3 bar (more moderate than 1.5)
        T_init=295,            # Further temp drop
        z_init=[0.5, 0.3, 0.2],
        F_init=10.0
    )
    
    # Add all streams to simulator
    for s in [s1, s2, s3, s4]:
        sim.add_stream(s)
    
    # =========================================================================
    # CREATE UNITS WITH IMPROVED SETTINGS
    # =========================================================================
    
    # FEED BOUNDARY - Source boundary condition
    feed = FeedBoundary(
        name='Feed',
        P_fixed=10e5,                # Fixed at 10 bar
        T_fixed=300,                 # Fixed at 300 K
        z_fixed=[0.5, 0.3, 0.2],    # Fixed composition
        compounds=compounds,
        F_fixed=10.0                 # Fixed flow rate 10 kg/s
    )
    
    # VALVE1 - Expansion valve (10 bar → 7 bar)
    # IMPROVED: Larger Cv and more moderate pressure drop
    valve1 = Valve(
        name='Valve1',
        Cv=12.0,                 # Larger Cv for smoother flow (was 8.0)
        opening=1.0,             # Fully open
        boundary_P=7e5,          # Target: 7 bar (was 8 bar)
        valve_type='globe'
    )
    
    # VOLUME - Pressure vessel/accumulator (2 m³)
    # IMPROVED: Larger volume for better stability
    volume = Volume(
        name='Volume_1',
        V_total=200.0,             # 2 m³ (was 1 m³, more stable)
        compounds=compounds,
        P_init=7e5,              # Initialize at 7 bar
        T_init=298,
        z_init=[0.5, 0.3, 0.2]
    )
    
    # VALVE2 - Relief valve (7 bar → 3 bar)
    # IMPROVED: Less extreme pressure drop
    valve2 = Valve(
        name='Valve2',
        Cv=8.0,                  # Reasonable Cv (was 5.0)
        opening=1.0,             # 90% open (was 80%, less restriction)
        boundary_P=3e5,          # Target: 3 bar (was 1.5 bar)
        valve_type='ball'
    )
    
    # PRODUCT BOUNDARY - Exit boundary condition
    product = Boundary(
        name='Product',
        P_fixed=3e5,                  # Fixed at 3 bar (was 1.5)
        T_fixed=300,                  # Fixed at 300 K (sink temperature)
        z_fixed=[0.5, 0.3, 0.2],     # Fixed composition
        compounds=compounds
    )
    
    # =========================================================================
    # CONNECT PORTS
    # =========================================================================
    
    # Feed → Valve1
    feed.connect_outlet(0, s1)
    valve1.connect_inlet(0, s1)
    
    # Valve1 → Volume
    valve1.connect_outlet(0, s2)
    volume.connect_inlet(0, s2)
    
    # Volume → Valve2
    volume.connect_outlet(0, s3)
    valve2.connect_inlet(0, s3)
    
    # Valve2 → Product
    valve2.connect_outlet(0, s4)
    product.connect_inlet(0, s4)
    
    # =========================================================================
    # ADD UNITS TO SIMULATOR
    # =========================================================================
    units = [feed, valve1, volume, valve2, product]
    
    for unit in units:
        # Load user code from master file
        unit.load_user_code_from_master_file()
        # Add to simulator
        sim.add_unit(unit)
    
    print("\n" + "="*80)
    print("IMPROVED FLOWSHEET BUILT SUCCESSFULLY")
    print("="*80)
    print(f"\nStreams:  {[s.name for s in [s1, s2, s3, s4]]}")
    print(f"Units:    {[u.name for u in units]}")
    print(f"\nSystem Architecture (IMPROVED):")
    print(f"  Feed (10 bar, 300K, 10 kg/s)")
    print(f"    ↓")
    print(f"  [Valve1: 10→7 bar, Cv=12.0, 100% open]")
    print(f"    ↓  (Expected ΔT: ~2K)")
    print(f"  [Volume: 2 m³ accumulator @ 7 bar]")
    print(f"    ↓")
    print(f"  [Valve2: 7→3 bar, Cv=8.0, 90% open]")
    print(f"    ↓  (Expected ΔT: ~3K)")
    print(f"  Product (3 bar, ~295K)")
    print(f"\nEXPECTED BEHAVIOR:")
    print(f"  • Total pressure drop: 10→3 bar (7 bar)")
    print(f"  • Total temp drop: ~5-6 K (300→294-295 K)")
    print(f"  • Flow maintained at ~10 kg/s")
    print(f"  • Volume stabilizes pressure at 7 bar")
    print(f"  • Less extreme than original (better stability)")
    print("="*80 + "\n")
    
    return sim


# ============================================================================
# ALTERNATIVE SCENARIOS TO TEST
# ============================================================================

def build_flowsheet_mild(sim, compounds):
    """
    SCENARIO 1: Very mild conditions (minimal pressure drop)
    Good for testing basic functionality
    
    10 bar → 9 bar → 9 bar → 8 bar → 8 bar
    Total drop: 2 bar, Expected ΔT: ~1-2K
    """
    s1 = Stream('s1', compounds, P_init=10e5, T_init=300, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s2 = Stream('s2', compounds, P_init=9e5, T_init=299.5, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s3 = Stream('s3', compounds, P_init=9e5, T_init=299.5, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s4 = Stream('s4', compounds, P_init=8e5, T_init=299, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    
    for s in [s1, s2, s3, s4]:
        sim.add_stream(s)
    
    feed = FeedBoundary('Feed', P_fixed=10e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], 
                        compounds=compounds, F_fixed=10.0)
    valve1 = Valve('Valve1', Cv=15.0, opening=1.0, boundary_P=9e5, valve_type='globe')
    volume = Volume('Volume_1', V_total=3.0, compounds=compounds, P_init=9e5, T_init=299.5, z_init=[0.5, 0.3, 0.2])
    valve2 = Valve('Valve2', Cv=10.0, opening=1.0, boundary_P=8e5, valve_type='ball')
    product = Boundary('Product', P_fixed=8e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], compounds=compounds)
    
    feed.connect_outlet(0, s1)
    valve1.connect_inlet(0, s1)
    valve1.connect_outlet(0, s2)
    volume.connect_inlet(0, s2)
    volume.connect_outlet(0, s3)
    valve2.connect_inlet(0, s3)
    valve2.connect_outlet(0, s4)
    product.connect_inlet(0, s4)
    
    units = [feed, valve1, volume, valve2, product]
    for unit in units:
        unit.load_user_code_from_master_file()
        sim.add_unit(unit)
    
    print("\nMILD SCENARIO: 10→9→8 bar (minimal drops)")
    return sim


def build_flowsheet_high_pressure(sim, compounds):
    """
    SCENARIO 2: High pressure system (all pressures elevated)
    Tests behavior at higher pressures where non-ideality matters
    
    50 bar → 40 bar → 40 bar → 30 bar → 30 bar
    Total drop: 20 bar, Expected ΔT: ~4-5K
    """
    s1 = Stream('s1', compounds, P_init=50e5, T_init=300, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s2 = Stream('s2', compounds, P_init=40e5, T_init=297, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s3 = Stream('s3', compounds, P_init=40e5, T_init=297, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s4 = Stream('s4', compounds, P_init=30e5, T_init=295, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    
    for s in [s1, s2, s3, s4]:
        sim.add_stream(s)
    
    feed = FeedBoundary('Feed', P_fixed=50e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], 
                        compounds=compounds, F_fixed=10.0)
    valve1 = Valve('Valve1', Cv=10.0, opening=1.0, boundary_P=40e5, valve_type='globe')
    volume = Volume('Volume_1', V_total=2.0, compounds=compounds, P_init=40e5, T_init=297, z_init=[0.5, 0.3, 0.2])
    valve2 = Valve('Valve2', Cv=8.0, opening=0.9, boundary_P=30e5, valve_type='ball')
    product = Boundary('Product', P_fixed=30e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], compounds=compounds)
    
    feed.connect_outlet(0, s1)
    valve1.connect_inlet(0, s1)
    valve1.connect_outlet(0, s2)
    volume.connect_inlet(0, s2)
    volume.connect_outlet(0, s3)
    valve2.connect_inlet(0, s3)
    valve2.connect_outlet(0, s4)
    product.connect_inlet(0, s4)
    
    units = [feed, valve1, volume, valve2, product]
    for unit in units:
        unit.load_user_code_from_master_file()
        sim.add_unit(unit)
    
    print("\nHIGH PRESSURE SCENARIO: 50→40→30 bar")
    return sim


def build_flowsheet_single_valve(sim, compounds):
    """
    SCENARIO 3: Simplified system (single valve, no volume)
    Good for debugging valve behavior in isolation
    
    Feed → Valve → Product
    10 bar → 5 bar
    """
    s1 = Stream('s1', compounds, P_init=10e5, T_init=300, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    s2 = Stream('s2', compounds, P_init=5e5, T_init=296, z_init=[0.5, 0.3, 0.2], F_init=10.0)
    
    for s in [s1, s2]:
        sim.add_stream(s)
    
    feed = FeedBoundary('Feed', P_fixed=10e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], 
                        compounds=compounds, F_fixed=10.0)
    valve = Valve('Valve1', Cv=10.0, opening=1.0, boundary_P=5e5, valve_type='globe')
    product = Boundary('Product', P_fixed=5e5, T_fixed=300, z_fixed=[0.5, 0.3, 0.2], compounds=compounds)
    
    feed.connect_outlet(0, s1)
    valve.connect_inlet(0, s1)
    valve.connect_outlet(0, s2)
    product.connect_inlet(0, s2)
    
    units = [feed, valve, product]
    for unit in units:
        unit.load_user_code_from_master_file()
        sim.add_unit(unit)
    
    print("\nSIMPLE SCENARIO: Feed→Valve→Product (10→5 bar)")
    return sim