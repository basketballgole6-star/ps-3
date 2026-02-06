"""
DAE (Differential-Algebraic Equation) Integrator
Handles systems with both differential and algebraic equations

System form:
    dy/dt = f(t, y, z)  - Differential equations
    0 = g(t, y, z)       - Algebraic equations

where y are differential variables, z are algebraic variables
"""

import numpy as np
from scipy.integrate import solve_ivp, ode
from typing import Callable, Tuple, Optional, Dict, List
import logging

logger = logging.getLogger(__name__)


class DAEIntegrator:
    """
    Integrator for Differential-Algebraic Equation systems
    
    Handles index-1 DAE systems commonly found in process simulation:
    - Differential equations: mass/energy balances in volumes
    - Algebraic equations: valve characteristics, stream properties
    """
    
    def __init__(self,
                 differential_equations: Callable,
                 algebraic_equations: Callable,
                 method: str = 'BDF',
                 rtol: float = 1e-6,
                 atol: float = 1e-8,
                 max_step: float = np.inf,
                 first_step: Optional[float] = None):
        """
        Args:
            differential_equations: Function f(t, y, z) returning dy/dt
            algebraic_equations: Function g(t, y, z) returning residual
            method: Integration method ('BDF', 'Radau', 'LSODA')
            rtol: Relative tolerance
            atol: Absolute tolerance
            max_step: Maximum time step
            first_step: Initial time step (None for automatic)
        """
        self.differential_equations = differential_equations
        self.algebraic_equations = algebraic_equations
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.first_step = first_step
        
        # Integration statistics
        self.n_steps = 0
        self.n_function_evals = 0
        self.n_jacobian_evals = 0
        self.n_algebraic_solves = 0
        
        # Storage for results
        self.t_history = []
        self.y_history = []
        self.z_history = []
        
    def integrate(self,
                  y0: np.ndarray,
                  z0: np.ndarray,
                  t_span: Tuple[float, float],
                  t_eval: Optional[np.ndarray] = None,
                  events: Optional[List[Callable]] = None) -> Dict:
        """
        Integrate DAE system from t0 to tf
        
        Args:
            y0: Initial values of differential variables
            z0: Initial values of algebraic variables
            t_span: (t0, tf) time span
            t_eval: Time points where solution is stored
            events: List of event functions for discontinuities
            
        Returns:
            Dictionary with solution arrays
        """
        self.t_history = []
        self.y_history = []
        self.z_history = []
        
        t0, tf = t_span
        y = y0.copy()
        z = z0.copy()
        
        logger.info(f"Starting DAE integration: t0={t0}, tf={tf}")
        logger.info(f"Differential vars: {len(y)}, Algebraic vars: {len(z)}")
        
        # Initialize solver based on method
        if self.method in ['BDF', 'Radau']:
            solution = self._integrate_scipy(y0, z0, t_span, t_eval, events)
        elif self.method == 'LSODA':
            solution = self._integrate_lsoda(y0, z0, t_span, t_eval)
        elif self.method == 'implicit_euler':
            solution = self._integrate_implicit_euler(y0, z0, t_span, t_eval)
        else:
            raise ValueError(f"Unknown integration method: {self.method}")
        
        logger.info(f"Integration completed: {self.n_steps} steps, "
                   f"{self.n_function_evals} function evaluations")
        
        return solution
    
    def _integrate_scipy(self,
                        y0: np.ndarray,
                        z0: np.ndarray,
                        t_span: Tuple[float, float],
                        t_eval: Optional[np.ndarray],
                        events: Optional[List[Callable]]) -> Dict:
        """
        Integrate using scipy.integrate.solve_ivp
        
        We reformulate the DAE as an ODE by solving algebraic equations at each step
        """
        t0, tf = t_span
        z_current = z0.copy()
        
        def ode_function(t, y):
            """
            Wrapper that solves algebraic equations, then returns dy/dt
            """
            nonlocal z_current
            
            # Solve algebraic equations: g(t, y, z) = 0
            from solver.newton_solver import NewtonSolver
            
            def algebraic_residual(z):
                return self.algebraic_equations(t, y, z)
            
            solver = NewtonSolver(max_iterations=20, tolerance=1e-6, verbose=False)
            z_new, converged, info = solver.solve(algebraic_residual, z_current)
            
            if not converged:
                logger.warning(f"Algebraic equations did not converge at t={t:.6f}")
            
            z_current = z_new
            self.n_algebraic_solves += 1
            
            # Compute dy/dt
            dydt = self.differential_equations(t, y, z_current)
            self.n_function_evals += 1
            
            return dydt
        
        # Event functions need access to z
        wrapped_events = None
        if events:
            wrapped_events = []
            for event_func in events:
                def wrapped_event(t, y, original_event=event_func):
                    return original_event(t, y, z_current)
                wrapped_events.append(wrapped_event)
        
        # Integrate
        sol = solve_ivp(
            ode_function,
            t_span,
            y0,
            method=self.method,
            t_eval=t_eval,
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
            first_step=self.first_step,
            events=wrapped_events,
            dense_output=True
        )
        
        self.n_steps = sol.nfev  # Approximate
        
        # Reconstruct z at each time point
        z_history = []
        for i in range(len(sol.t)):
            t = sol.t[i]
            y = sol.y[:, i]
            
            # Solve for z at this point
            def algebraic_residual(z):
                return self.algebraic_equations(t, y, z)
            
            from solver.newton_solver import NewtonSolver
            solver = NewtonSolver(max_iterations=20, tolerance=1e-6, verbose=False)
            z_sol, _, _ = solver.solve(algebraic_residual, z_current)
            z_current = z_sol
            z_history.append(z_sol)
        
        return {
            't': sol.t,
            'y': sol.y,
            'z': np.array(z_history).T,
            'success': sol.success,
            'message': sol.message,
            'nfev': self.n_function_evals,
            'njev': self.n_jacobian_evals,
        }
    
    def _integrate_implicit_euler(self,
                                   y0: np.ndarray,
                                   z0: np.ndarray,
                                   t_span: Tuple[float, float],
                                   t_eval: Optional[np.ndarray]) -> Dict:
        """
        Implicit Euler method for DAE (more robust for stiff systems)
        
        Solves at each step:
            y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1}, z_{n+1})
            0 = g(t_{n+1}, y_{n+1}, z_{n+1})
        """
        from solver.newton_solver import NewtonSolver
        
        t0, tf = t_span
        
        if t_eval is None:
            # Automatic time stepping
            dt_initial = self.first_step if self.first_step else (tf - t0) / 100
            t_eval = self._generate_time_grid(t0, tf, dt_initial)
        
        n_points = len(t_eval)
        n_y = len(y0)
        n_z = len(z0)
        
        # Allocate storage
        y_solution = np.zeros((n_y, n_points))
        z_solution = np.zeros((n_z, n_points))
        
        y_solution[:, 0] = y0
        z_solution[:, 0] = z0
        
        y_current = y0.copy()
        z_current = z0.copy()
        
        for i in range(1, n_points):
            t_old = t_eval[i-1]
            t_new = t_eval[i]
            dt = t_new - t_old
            
            # Solve coupled system for y_new and z_new
            def residual(x):
                y_new = x[:n_y]
                z_new = x[n_y:]
                
                # Differential equation residual
                dydt = self.differential_equations(t_new, y_new, z_new)
                residual_diff = y_new - y_current - dt * dydt
                
                # Algebraic equation residual
                residual_alg = self.algebraic_equations(t_new, y_new, z_new)
                
                return np.concatenate([residual_diff, residual_alg])
            
            # Initial guess
            x0 = np.concatenate([y_current, z_current])
            
            # Solve
            solver = NewtonSolver(max_iterations=50, tolerance=1e-6, verbose=False)
            x_solution, converged, info = solver.solve(residual, x0)
            
            if not converged:
                logger.warning(f"Implicit Euler did not converge at t={t_new:.6f}")
            
            y_current = x_solution[:n_y]
            z_current = x_solution[n_y:]
            
            y_solution[:, i] = y_current
            z_solution[:, i] = z_current
            
            self.n_steps += 1
            self.n_function_evals += info['iterations'] * 2  # Approximate
        
        return {
            't': t_eval,
            'y': y_solution,
            'z': z_solution,
            'success': True,
            'message': 'Integration completed',
            'nfev': self.n_function_evals,
            'njev': 0,
        }
    
    def _integrate_lsoda(self,
                        y0: np.ndarray,
                        z0: np.ndarray,
                        t_span: Tuple[float, float],
                        t_eval: Optional[np.ndarray]) -> Dict:
        """
        LSODA integration (switches between stiff/non-stiff automatically)
        Similar to scipy method but with LSODA specifically
        """
        # Similar implementation to _integrate_scipy but using ode class
        # with 'lsoda' integrator
        
        t0, tf = t_span
        z_current = z0.copy()
        
        def ode_function(t, y):
            nonlocal z_current
            
            # Solve algebraic equations
            from solver.newton_solver import NewtonSolver
            
            def algebraic_residual(z):
                return self.algebraic_equations(t, y, z)
            
            solver = NewtonSolver(max_iterations=20, tolerance=1e-6, verbose=False)
            z_new, _, _ = solver.solve(algebraic_residual, z_current)
            z_current = z_new
            
            # Return dy/dt
            return self.differential_equations(t, y, z_current)
        
        integrator = ode(ode_function)
        integrator.set_integrator('lsoda', rtol=self.rtol, atol=self.atol, max_step=self.max_step)
        integrator.set_initial_value(y0, t0)
        
        if t_eval is None:
            t_eval = np.linspace(t0, tf, 100)
        
        y_solution = np.zeros((len(y0), len(t_eval)))
        z_solution = []
        
        y_solution[:, 0] = y0
        z_solution.append(z0)
        
        for i, t in enumerate(t_eval[1:], 1):
            integrator.integrate(t)
            y_solution[:, i] = integrator.y
            z_solution.append(z_current.copy())
        
        return {
            't': t_eval,
            'y': y_solution,
            'z': np.array(z_solution).T,
            'success': integrator.successful(),
            'message': 'Integration completed',
            'nfev': self.n_function_evals,
            'njev': 0,
        }
    
    def _generate_time_grid(self, t0: float, tf: float, dt_initial: float) -> np.ndarray:
        """Generate adaptive time grid"""
        n_steps = int(np.ceil((tf - t0) / dt_initial))
        return np.linspace(t0, tf, n_steps + 1)
    
    def get_statistics(self) -> Dict:
        """Return integration statistics"""
        return {
            'n_steps': self.n_steps,
            'n_function_evals': self.n_function_evals,
            'n_jacobian_evals': self.n_jacobian_evals,
            'n_algebraic_solves': self.n_algebraic_solves,
        }


class AdaptiveDAEIntegrator(DAEIntegrator):
    """
    DAE integrator with adaptive step size control
    """
    
    def __init__(self, *args, min_dt: float = 1e-6, max_dt: float = 1.0,
                 safety_factor: float = 0.9, **kwargs):
        super().__init__(*args, **kwargs)
        self.min_dt = min_dt
        self.max_dt = max_dt
        self.safety_factor = safety_factor
    
    def integrate(self, y0, z0, t_span, t_eval=None, events=None):
        """
        Integrate with adaptive time stepping
        """
        if self.method == 'adaptive_euler':
            return self._integrate_adaptive_euler(y0, z0, t_span, t_eval)
        else:
            # Use parent class method
            return super().integrate(y0, z0, t_span, t_eval, events)
    
    def _integrate_adaptive_euler(self, y0, z0, t_span, t_eval):
        """
        Adaptive Implicit Euler with error estimation
        """
        from solver.newton_solver import NewtonSolver
        
        t0, tf = t_span
        
        # Initial conditions
        t = t0
        y = y0.copy()
        z = z0.copy()
        dt = self.first_step if self.first_step else (tf - t0) / 100
        
        # Storage
        t_history = [t0]
        y_history = [y0]
        z_history = [z0]
        
        while t < tf:
            # Don't overshoot
            dt = min(dt, tf - t)
            
            # Take step with dt
            y1, z1, error1 = self._implicit_euler_step(t, y, z, dt)
            
            # Take two steps with dt/2
            y_half, z_half, _ = self._implicit_euler_step(t, y, z, dt/2)
            y2, z2, error2 = self._implicit_euler_step(t + dt/2, y_half, z_half, dt/2)
            
            # Estimate error (Richardson extrapolation)
            error_estimate = np.linalg.norm(y2 - y1) / np.linalg.norm(y2 + 1e-10)
            
            # Accept or reject step
            if error_estimate < self.rtol or dt <= self.min_dt:
                # Accept step
                t += dt
                y = y2  # Use more accurate solution
                z = z2
                
                t_history.append(t)
                y_history.append(y)
                z_history.append(z)
                
                self.n_steps += 1
            
            # Adjust step size
            if error_estimate > 0:
                dt_new = self.safety_factor * dt * (self.rtol / error_estimate) ** 0.5
                dt = np.clip(dt_new, self.min_dt, self.max_dt)
            
            if t >= tf:
                break
        
        # Interpolate to t_eval if provided
        if t_eval is not None:
            from scipy.interpolate import interp1d
            y_array = np.array(y_history).T
            z_array = np.array(z_history).T
            
            y_interp = interp1d(t_history, y_array, axis=1, kind='cubic')
            z_interp = interp1d(t_history, z_array, axis=1, kind='cubic')
            
            return {
                't': t_eval,
                'y': y_interp(t_eval),
                'z': z_interp(t_eval),
                'success': True,
                'message': 'Adaptive integration completed',
            }
        else:
            return {
                't': np.array(t_history),
                'y': np.array(y_history).T,
                'z': np.array(z_history).T,
                'success': True,
                'message': 'Adaptive integration completed',
            }
    
    def _implicit_euler_step(self, t, y, z, dt):
        """
        Single implicit Euler step with error estimate
        """
        from solver.newton_solver import NewtonSolver
        
        n_y = len(y)
        n_z = len(z)
        
        def residual(x):
            y_new = x[:n_y]
            z_new = x[n_y:]
            
            dydt = self.differential_equations(t + dt, y_new, z_new)
            residual_diff = y_new - y - dt * dydt
            
            residual_alg = self.algebraic_equations(t + dt, y_new, z_new)
            
            return np.concatenate([residual_diff, residual_alg])
        
        x0 = np.concatenate([y, z])
        
        solver = NewtonSolver(max_iterations=50, tolerance=1e-8, verbose=False)
        x_solution, converged, info = solver.solve(residual, x0)
        
        y_new = x_solution[:n_y]
        z_new = x_solution[n_y:]
        error = info['final_residual']
        
        return y_new, z_new, error


if __name__ == "__main__":
    # Test: Simple DAE system
    # Differential: dy/dt = -y + z
    # Algebraic: 0 = y - z^2
    # Exact solution: y(t) = exp(-t), z(t) = sqrt(exp(-t))
    
    def diff_eqs(t, y, z):
        return np.array([-y[0] + z[0]])
    
    def alg_eqs(t, y, z):
        return np.array([y[0] - z[0]**2])
    
    y0 = np.array([1.0])
    z0 = np.array([1.0])
    t_span = (0, 5)
    t_eval = np.linspace(0, 5, 51)
    
    integrator = DAEIntegrator(diff_eqs, alg_eqs, method='BDF')
    solution = integrator.integrate(y0, z0, t_span, t_eval)
    
    # Plot
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(10, 4))
    
    plt.subplot(1, 2, 1)
    plt.plot(solution['t'], solution['y'][0], 'o-', label='y(t) - computed')
    plt.plot(solution['t'], np.exp(-solution['t']), '--', label='y(t) - exact')
    plt.xlabel('Time')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.plot(solution['t'], solution['z'][0], 'o-', label='z(t) - computed')
    plt.plot(solution['t'], np.sqrt(np.exp(-solution['t'])), '--', label='z(t) - exact')
    plt.xlabel('Time')
    plt.ylabel('z')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nIntegration successful: {solution['success']}")
    print(f"Message: {solution['message']}")
    print(f"Function evaluations: {solution['nfev']}")