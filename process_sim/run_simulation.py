# run_simulation.py

from config.components import COMPONENTS
from solver.pf_solver import DynamicSimulator
from flowsheet.flowsheet import build_flowsheet

from plot_results import ask_and_plot


import config.sim_settings as sim_settings

def main():
    # Simulation parameters
    t_total = 0.02   # total simulation time [s] 
    dt = 0.02         # time step [s] 
    dt_result = 0.02   # result storage interval [s] 
 
    # --- Add these lines to update global settings --- 
    sim_settings.TOTAL_TIME = t_total 
    sim_settings.DT = dt 
    sim_settings.RT = dt_result 
 
    # Create simulator and build flowsheet 
    sim = DynamicSimulator(COMPONENTS) 
    sim = build_flowsheet(sim, COMPONENTS) 
 
    # Run simulation 
    time_points, iterations_history, convergence_flags = sim.run(t_total, dt, dt_result) 

    # Plot results

    # Print final results summary
    print("\n" + "="*70)
    print("FINAL RESULTS")
    print("="*70)
    for volume in sim.volumes:
        F_in = sum((port.stream.F for port in volume.inlet_ports if port.stream is not None))
        F_out = sum((port.stream.F for port in volume.outlet_ports if port.stream is not None))
        print(f"{volume.name}: P={volume.P/1e5:.2f} bar, F_in={F_in:.2f}, F_out={F_out:.2f}, Balance={F_in-F_out:.4f}")
    print("="*70 + "\n")
    ask_and_plot(sim, time_points)

if __name__ == "__main__":
    main()
