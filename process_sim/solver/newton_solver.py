"""
Newton-Raphson Solver for Coupled Nonlinear Equations
Solves systems of the form F(x) = 0 using Newton's method
"""

import numpy as np
import logging
from typing import Callable, Tuple, Optional, List, Dict

logger = logging.getLogger(__name__)


class NewtonSolver:
    """
    Newton-Raphson solver for nonlinear equation systems.
    
    Solves: F(x) = 0
    Using: x_{n+1} = x_n - J^{-1} * F(x_n)
    where J is the Jacobian matrix
    """
    
    def __init__(self, 
                 max_iterations: int = 50,
                 tolerance: float = 1e-6,
                 min_step: float = 1e-10,
                 max_step: float = 1.0,
                 damping: bool = True,
                 verbose: bool = False):
        """
        Args:
            max_iterations: Maximum number of Newton iterations
            tolerance: Convergence tolerance (||F(x)|| < tolerance)
            min_step: Minimum allowable step size
            max_step: Maximum allowable step size (for damping)
            damping: Use line search damping
            verbose: Print iteration details
        """
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.min_step = min_step
        self.max_step = max_step
        self.damping = damping
        self.verbose = verbose
        
        # Diagnostics
        self.iteration_history = []
        self.residual_history = []
        self.step_history = []
        
    def solve(self, 
              residual_function: Callable,
              x0: np.ndarray,
              jacobian_function: Optional[Callable] = None,
              bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> Tuple[np.ndarray, bool, Dict]:
        """
        Solve nonlinear system F(x) = 0
        
        Args:
            residual_function: Function that returns F(x) given x
            x0: Initial guess
            jacobian_function: Function that returns Jacobian J(x). If None, uses finite differences
            bounds: Optional (lower, upper) bounds for variables
            
        Returns:
            solution: Final x value
            converged: True if converged
            info: Dictionary with convergence information
        """
        x = np.array(x0, dtype=float)
        n_vars = len(x)
        
        # Reset diagnostics
        self.iteration_history = []
        self.residual_history = []
        self.step_history = []
        
        if self.verbose:
            logger.info(f"Starting Newton solver: {n_vars} variables")
        
        # Initial residual
        F = residual_function(x)
        residual_norm = np.linalg.norm(F)
        self.residual_history.append(residual_norm)
        
        if self.verbose:
            logger.info(f"Initial residual: {residual_norm:.6e}")
        
        # Check if already converged
        if residual_norm < self.tolerance:
            return x, True, self._build_info(0, True, "Already converged")
        
        # Newton iterations
        for iteration in range(self.max_iterations):
            
            # Compute Jacobian
            if jacobian_function is not None:
                J = jacobian_function(x)
            else:
                J = self._finite_difference_jacobian(residual_function, x, F)
            
            # Check Jacobian conditioning
            try:
                cond = np.linalg.cond(J)
                if cond > 1e12:
                    logger.warning(f"Iteration {iteration}: Jacobian ill-conditioned (cond={cond:.2e})")
            except:
                pass
            
            # Solve linear system: J * dx = -F
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                logger.error(f"Iteration {iteration}: Singular Jacobian")
                return x, False, self._build_info(iteration, False, "Singular Jacobian")
            
            # Apply step limiting
            step_norm = np.linalg.norm(dx)
            if step_norm > self.max_step:
                dx = dx * (self.max_step / step_norm)
                step_norm = self.max_step
            
            # Line search with damping
            if self.damping:
                x_new, F_new, alpha = self._line_search(residual_function, x, F, dx, residual_norm)
            else:
                alpha = 1.0
                x_new = x + dx
                if bounds is not None:
                    x_new = self._apply_bounds(x_new, bounds)
                F_new = residual_function(x_new)
            
            # Update
            x = x_new
            F = F_new
            residual_norm_new = np.linalg.norm(F)
            
            # Store history
            self.iteration_history.append(iteration)
            self.residual_history.append(residual_norm_new)
            self.step_history.append(step_norm * alpha)
            
            if self.verbose:
                logger.info(f"Iter {iteration:3d}: residual={residual_norm_new:.6e}, "
                          f"step={step_norm*alpha:.6e}, alpha={alpha:.3f}")
            
            # Check convergence
            if residual_norm_new < self.tolerance:
                if self.verbose:
                    logger.info(f"Converged in {iteration+1} iterations")
                return x, True, self._build_info(iteration+1, True, "Converged")
            
            # Check for stagnation
            if step_norm * alpha < self.min_step:
                logger.warning(f"Iteration {iteration}: Step size too small, possible stagnation")
                return x, False, self._build_info(iteration+1, False, "Step size too small")
            
            # Check if residual is increasing significantly
            if residual_norm_new > 10 * residual_norm:
                logger.warning(f"Iteration {iteration}: Residual increasing significantly")
                # Continue anyway, might recover
            
            residual_norm = residual_norm_new
        
        # Did not converge
        logger.warning(f"Newton solver did not converge in {self.max_iterations} iterations")
        return x, False, self._build_info(self.max_iterations, False, "Max iterations reached")
    
    def _finite_difference_jacobian(self, 
                                    residual_function: Callable,
                                    x: np.ndarray,
                                    F: np.ndarray,
                                    epsilon: float = 1e-8) -> np.ndarray:
        """
        Compute Jacobian using forward finite differences
        J[i,j] = dF_i/dx_j ≈ (F_i(x + h*e_j) - F_i(x)) / h
        """
        n = len(x)
        J = np.zeros((n, n))
        
        for j in range(n):
            # Perturb x[j]
            x_perturbed = x.copy()
            h = epsilon * max(abs(x[j]), 1.0)
            x_perturbed[j] += h
            
            # Compute perturbed residual
            F_perturbed = residual_function(x_perturbed)
            
            # Finite difference
            J[:, j] = (F_perturbed - F) / h
        
        return J
    
    def _line_search(self,
                    residual_function: Callable,
                    x: np.ndarray,
                    F: np.ndarray,
                    dx: np.ndarray,
                    residual_norm: float,
                    max_backtracks: int = 10) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Backtracking line search to ensure residual decreases
        
        Returns:
            x_new: New point
            F_new: Residual at new point
            alpha: Step size used
        """
        alpha = 1.0
        rho = 0.5  # Backtracking factor
        c = 1e-4   # Armijo constant
        
        # Target reduction (Armijo condition)
        target = residual_norm**2 * (1 - c * alpha)
        
        for backtrack in range(max_backtracks):
            x_new = x + alpha * dx
            F_new = residual_function(x_new)
            residual_norm_new = np.linalg.norm(F_new)
            
            # Check if sufficient decrease
            if residual_norm_new**2 <= target:
                return x_new, F_new, alpha
            
            # Reduce step size
            alpha *= rho
            target = residual_norm**2 * (1 - c * alpha)
        
        # If line search fails, return original step
        logger.warning("Line search failed, using full step")
        x_new = x + dx
        F_new = residual_function(x_new)
        return x_new, F_new, 1.0
    
    def _apply_bounds(self, x: np.ndarray, bounds: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
        """Apply bounds to variables"""
        lower, upper = bounds
        return np.clip(x, lower, upper)
    
    def _build_info(self, iterations: int, converged: bool, message: str) -> Dict:
        """Build information dictionary"""
        return {
            'iterations': iterations,
            'converged': converged,
            'message': message,
            'final_residual': self.residual_history[-1] if self.residual_history else None,
            'residual_history': self.residual_history.copy(),
            'step_history': self.step_history.copy(),
        }
    
    def plot_convergence(self):
        """Plot convergence history"""
        import matplotlib.pyplot as plt
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
        # Residual history
        ax1.semilogy(self.iteration_history, self.residual_history, 'o-')
        ax1.axhline(self.tolerance, color='r', linestyle='--', label='Tolerance')
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('Residual Norm')
        ax1.set_title('Convergence History')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Step size history
        ax2.semilogy(self.iteration_history, self.step_history, 'o-')
        ax2.axhline(self.min_step, color='r', linestyle='--', label='Min step')
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('Step Size')
        ax2.set_title('Step Size History')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.show()


class ModifiedNewtonSolver(NewtonSolver):
    """
    Modified Newton solver with additional features:
    - Trust region method
    - Adaptive damping
    - Better handling of ill-conditioned systems
    """
    
    def __init__(self, *args, trust_region: float = 1.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.trust_region = trust_region
    
    def solve(self, residual_function, x0, jacobian_function=None, bounds=None):
        """Enhanced Newton solver with trust region"""
        x = np.array(x0, dtype=float)
        n_vars = len(x)
        
        self.iteration_history = []
        self.residual_history = []
        self.step_history = []
        
        F = residual_function(x)
        residual_norm = np.linalg.norm(F)
        self.residual_history.append(residual_norm)
        
        if residual_norm < self.tolerance:
            return x, True, self._build_info(0, True, "Already converged")
        
        trust_radius = self.trust_region
        
        for iteration in range(self.max_iterations):
            
            # Compute Jacobian
            if jacobian_function is not None:
                J = jacobian_function(x)
            else:
                J = self._finite_difference_jacobian(residual_function, x, F)
            
            # Solve trust region subproblem
            # min ||F + J*dx||^2  subject to ||dx|| <= trust_radius
            
            # Try full Newton step
            try:
                dx_newton = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                # Singular Jacobian - use damped least squares
                lambda_damp = 1e-6
                J_damped = J.T @ J + lambda_damp * np.eye(n_vars)
                dx_newton = np.linalg.solve(J_damped, -J.T @ F)
            
            step_norm = np.linalg.norm(dx_newton)
            
            # Apply trust region
            if step_norm <= trust_radius:
                dx = dx_newton
            else:
                # Scale to trust region boundary
                dx = dx_newton * (trust_radius / step_norm)
            
            # Trial step
            x_new = x + dx
            if bounds is not None:
                x_new = self._apply_bounds(x_new, bounds)
            
            F_new = residual_function(x_new)
            residual_norm_new = np.linalg.norm(F_new)
            
            # Compute reduction ratio
            # rho = (actual reduction) / (predicted reduction)
            predicted_reduction = residual_norm**2 - np.linalg.norm(F + J @ dx)**2
            actual_reduction = residual_norm**2 - residual_norm_new**2
            
            if predicted_reduction > 0:
                rho = actual_reduction / predicted_reduction
            else:
                rho = 0
            
            # Update trust region
            if rho < 0.25:
                # Bad step - shrink trust region
                trust_radius *= 0.25
            elif rho > 0.75 and step_norm >= 0.9 * trust_radius:
                # Good step at boundary - expand trust region
                trust_radius = min(2 * trust_radius, self.max_step)
            
            # Accept or reject step
            if rho > 0.1:  # Accept
                x = x_new
                F = F_new
                residual_norm = residual_norm_new
            else:  # Reject
                logger.debug(f"Iteration {iteration}: Step rejected (rho={rho:.3f})")
            
            # Store history
            self.iteration_history.append(iteration)
            self.residual_history.append(residual_norm)
            self.step_history.append(np.linalg.norm(dx))
            
            if self.verbose:
                logger.info(f"Iter {iteration:3d}: residual={residual_norm:.6e}, "
                          f"step={np.linalg.norm(dx):.6e}, trust={trust_radius:.3e}, rho={rho:.3f}")
            
            # Check convergence
            if residual_norm < self.tolerance:
                if self.verbose:
                    logger.info(f"Converged in {iteration+1} iterations")
                return x, True, self._build_info(iteration+1, True, "Converged")
            
            # Check for stagnation
            if trust_radius < self.min_step:
                logger.warning("Trust region too small, possible stagnation")
                return x, False, self._build_info(iteration+1, False, "Trust region too small")
        
        logger.warning(f"Modified Newton solver did not converge in {self.max_iterations} iterations")
        return x, False, self._build_info(self.max_iterations, False, "Max iterations reached")


# Convenience function
def solve_nonlinear_system(residual_function: Callable,
                          x0: np.ndarray,
                          method: str = 'newton',
                          **kwargs) -> Tuple[np.ndarray, bool, Dict]:
    """
    Convenience function to solve nonlinear system
    
    Args:
        residual_function: Function F(x) to find root of
        x0: Initial guess
        method: 'newton' or 'modified_newton'
        **kwargs: Additional arguments for solver
        
    Returns:
        solution, converged, info
    """
    if method == 'newton':
        solver = NewtonSolver(**kwargs)
    elif method == 'modified_newton':
        solver = ModifiedNewtonSolver(**kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return solver.solve(residual_function, x0)


if __name__ == "__main__":
    # Test: Solve simple system
    # F1 = x^2 + y^2 - 4 = 0
    # F2 = x*y - 1 = 0
    # Solution: (x,y) ≈ (1.732, 0.577) or (-1.732, -0.577)
    
    def residual(x):
        return np.array([
            x[0]**2 + x[1]**2 - 4,
            x[0] * x[1] - 1
        ])
    
    x0 = np.array([2.0, 1.0])
    
    solver = NewtonSolver(verbose=True)
    solution, converged, info = solver.solve(residual, x0)
    
    print(f"\nSolution: {solution}")
    print(f"Converged: {converged}")
    print(f"Iterations: {info['iterations']}")
    print(f"Final residual: {info['final_residual']:.6e}")
    print(f"Verification: F(x) = {residual(solution)}")
    
    solver.plot_convergence()