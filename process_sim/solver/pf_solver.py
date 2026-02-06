# solver/pf_solver.py
import config.sim_settings as sim_settings
class DynamicSimulator:
    """
    Manages the process simulation: time stepping, convergence, and results.
    """
    def __init__(self, compounds):
        self.compounds = compounds
        self.streams = {}
        self.units = []
        self.volumes = []  # For history/plotting

    def add_stream(self, stream):
        self.streams[stream.name] = stream

    def add_unit(self, unit):
        self.units.append(unit)
        # If it's a Volume, keep for plotting
        if hasattr(unit, 'P_history'):
            self.volumes.append(unit)

    def solve_timestep(self, dt, current_time=None):
        max_iter = 100
        tolerance = 1e-3
        relaxation_P = 0.4
        relaxation_F = 0.3

        for iteration in range(max_iter):
            old_values = {}
            for name, stream in self.streams.items():
                old_values[name] = (stream.P, stream.T, stream.F)

            # 1. Calculate all units except volumes
            for unit in self.units:
                # Volumes need dt, others don't
                if hasattr(unit, 'calculate') and not hasattr(unit, 'P_history'):
                    unit.calculate()
            # 2. Calculate all volumes (with dt)
            for unit in self.units:
                if hasattr(unit, 'P_history'):
                    unit.calculate(dt)

            # Relaxation
            if iteration > 0:
                for name, stream in self.streams.items():
                    old_P, old_T, old_F = old_values[name]
                    stream.P = relaxation_P * stream.P + (1 - relaxation_P) * old_P
                    stream.T = relaxation_P * stream.T + (1 - relaxation_P) * old_T
                    stream.F = relaxation_F * stream.F + (1 - relaxation_F) * old_F

            # Check convergence
            max_change = 0
            for name, stream in self.streams.items():
                old_P, old_T, old_F = old_values[name]
                change_P = abs(stream.P - old_P) / (abs(old_P) + 1e3)
                change_T = abs(stream.T - old_T) / (abs(old_T) + 1.0)
                change_F = abs(stream.F - old_F) / (abs(old_F) + 1e-3)
                max_change = max(max_change, change_P, change_T, change_F)

            if max_change < tolerance:
                return iteration + 1, True

        return max_iter, False

    def run(self, t_total, dt, dt_result):
        n_steps = int(t_total / dt)
        n_result_steps = int(dt_result / dt)

        time_points = []
        iterations_history = []
        convergence_flags = []

        print(f"\n{'='*70}")
        print(f"DYNAMIC PROCESS SIMULATOR")
        print(f"{'='*70}")
        print(f"Total time: {t_total} s")
        print(f"Time step: {dt} s")
        print(f"Result interval: {dt_result} s")
        print(f"Total steps: {n_steps}")
        print(f"{'='*70}\n")

        time_points.append(0.0)
        for stream in self.streams.values():
            stream.store_history()
        for volume in self.volumes:
            if hasattr(volume, 'store_history'):
                volume.store_history()

        converged_count = 0

        for step in range(n_steps):
            t = (step + 1) * dt
            sim_settings.CT = t  # Update global CT
            print("CT in simulator loop:", sim_settings.CT)  # <--- Add this line
            n_iter, converged = self.solve_timestep(dt, current_time=t)
            iterations_history.append(n_iter)
            convergence_flags.append(converged)

            if converged:
                converged_count += 1

            if (step + 1) % n_result_steps == 0:
                time_points.append(t)
                for stream in self.streams.values():
                    stream.run_user_code()
                    stream.store_history()
                for unit in self.units:
                    unit.run_user_code('before')  # or 'after', as appropriate
                    unit.store_history()
                for volume in self.volumes:
                    if hasattr(volume, 'store_history'):
                        volume.store_history()
                status = "✓" if converged else "✗"
               
                for volume in self.volumes:
                    F_in = sum((port.stream.F for port in volume.inlet_ports if port.stream is not None))
                    F_out = sum((port.stream.F for port in volume.outlet_ports if port.stream is not None))
                

        convergence_rate = 100 * converged_count / n_steps

        print(f"\n{'='*70}")
        print(f"Simulation Complete!")
        print(f"Convergence rate: {convergence_rate:.1f}%")
        print(f"{'='*70}\n")

        return time_points, iterations_history, convergence_flags
