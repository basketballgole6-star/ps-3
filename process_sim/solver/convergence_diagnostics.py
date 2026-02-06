"""
Convergence Diagnostics and Analysis Tools
Helps debug why simulations don't converge
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class ConvergenceDiagnostics:
    """
    Analyze and diagnose convergence issues in process simulations
    """
    
    def __init__(self):
        self.iteration_data = []
        self.variable_history = {}
        self.residual_history = {}
        self.equation_names = []
        self.variable_names = []
        
    def record_iteration(self,
                        iteration: int,
                        variables: np.ndarray,
                        residuals: np.ndarray,
                        variable_names: Optional[List[str]] = None,
                        equation_names: Optional[List[str]] = None):
        """
        Record data from an iteration
        
        Args:
            iteration: Iteration number
            variables: Current variable values
            residuals: Current residual values
            variable_names: Names of variables
            equation_names: Names of equations
        """
        if variable_names and not self.variable_names:
            self.variable_names = variable_names
        if equation_names and not self.equation_names:
            self.equation_names = equation_names
        
        # Store iteration data
        self.iteration_data.append({
            'iteration': iteration,
            'variables': variables.copy(),
            'residuals': residuals.copy(),
            'residual_norm': np.linalg.norm(residuals),
        })
        
        # Store variable history
        for i, var in enumerate(variables):
            var_name = self.variable_names[i] if i < len(self.variable_names) else f"var_{i}"
            if var_name not in self.variable_history:
                self.variable_history[var_name] = []
            self.variable_history[var_name].append((iteration, var))
        
        # Store residual history
        for i, res in enumerate(residuals):
            eq_name = self.equation_names[i] if i < len(self.equation_names) else f"eq_{i}"
            if eq_name not in self.residual_history:
                self.residual_history[eq_name] = []
            self.residual_history[eq_name].append((iteration, res))
    
    def analyze_convergence(self, tolerance: float = 1e-6) -> Dict:
        """
        Analyze convergence behavior and identify issues
        
        Returns:
            Dictionary with diagnostic information
        """
        if not self.iteration_data:
            return {'error': 'No iteration data recorded'}
        
        diagnostics = {}
        
        # 1. Overall convergence status
        final_residual = self.iteration_data[-1]['residual_norm']
        converged = final_residual < tolerance
        diagnostics['converged'] = converged
        diagnostics['final_residual'] = final_residual
        diagnostics['n_iterations'] = len(self.iteration_data)
        
        # 2. Convergence rate
        if len(self.iteration_data) > 2:
            residual_norms = [data['residual_norm'] for data in self.iteration_data]
            rate = self._estimate_convergence_rate(residual_norms)
            diagnostics['convergence_rate'] = rate
            
            if rate > 0.9:
                diagnostics['warning'] = 'Slow convergence (rate > 0.9)'
            elif rate > 1.0:
                diagnostics['warning'] = 'Diverging (rate > 1.0)'
        
        # 3. Identify problematic equations
        problematic_equations = self._identify_problematic_equations(tolerance)
        if problematic_equations:
            diagnostics['problematic_equations'] = problematic_equations
        
        # 4. Identify oscillating variables
        oscillating_vars = self._identify_oscillating_variables()
        if oscillating_vars:
            diagnostics['oscillating_variables'] = oscillating_vars
        
        # 5. Identify stagnant variables
        stagnant_vars = self._identify_stagnant_variables()
        if stagnant_vars:
            diagnostics['stagnant_variables'] = stagnant_vars
        
        # 6. Check for ill-conditioning
        if len(self.iteration_data) > 1:
            last_vars = self.iteration_data[-1]['variables']
            prev_vars = self.iteration_data[-2]['variables']
            var_change = np.abs(last_vars - prev_vars)
            res_change = np.abs(self.iteration_data[-1]['residuals'] - 
                               self.iteration_data[-2]['residuals'])
            
            # If small variable change causes large residual change → ill-conditioned
            if np.any(var_change < 1e-6) and np.any(res_change > 1e-3):
                diagnostics['warning_ill_conditioned'] = True
        
        # 7. Suggest improvements
        suggestions = self._generate_suggestions(diagnostics)
        diagnostics['suggestions'] = suggestions
        
        return diagnostics
    
    def _estimate_convergence_rate(self, residual_norms: List[float]) -> float:
        """
        Estimate convergence rate from residual history
        Linear convergence: ||r_{n+1}|| ≈ C * ||r_n||
        Rate C should be < 1 for convergence
        """
        if len(residual_norms) < 3:
            return None
        
        # Use last few iterations
        recent = residual_norms[-5:]
        rates = []
        
        for i in range(1, len(recent)):
            if recent[i-1] > 1e-15:  # Avoid division by zero
                rate = recent[i] / recent[i-1]
                rates.append(rate)
        
        return np.mean(rates) if rates else None
    
    def _identify_problematic_equations(self, tolerance: float) -> List[Tuple[str, float]]:
        """
        Identify equations with large residuals
        """
        if not self.iteration_data:
            return []
        
        final_residuals = self.iteration_data[-1]['residuals']
        problematic = []
        
        for i, res in enumerate(final_residuals):
            if abs(res) > tolerance * 10:  # 10x tolerance
                eq_name = self.equation_names[i] if i < len(self.equation_names) else f"eq_{i}"
                problematic.append((eq_name, res))
        
        # Sort by magnitude
        problematic.sort(key=lambda x: abs(x[1]), reverse=True)
        
        return problematic
    
    def _identify_oscillating_variables(self, threshold: float = 0.1) -> List[str]:
        """
        Identify variables that are oscillating between iterations
        """
        oscillating = []
        
        for var_name, history in self.variable_history.items():
            if len(history) < 4:
                continue
            
            # Get recent values
            recent_values = [val for _, val in history[-6:]]
            
            # Check for oscillation: sign of change alternates
            changes = np.diff(recent_values)
            sign_changes = np.diff(np.sign(changes))
            
            # If many sign changes → oscillating
            n_sign_changes = np.sum(np.abs(sign_changes) > 0)
            if n_sign_changes >= len(changes) - 1:
                # Also check magnitude
                avg_change = np.mean(np.abs(changes))
                avg_value = np.mean(np.abs(recent_values))
                if avg_change > threshold * avg_value:
                    oscillating.append(var_name)
        
        return oscillating
    
    def _identify_stagnant_variables(self, threshold: float = 1e-10) -> List[str]:
        """
        Identify variables that are not changing
        """
        stagnant = []
        
        for var_name, history in self.variable_history.items():
            if len(history) < 3:
                continue
            
            # Get recent values
            recent_values = [val for _, val in history[-5:]]
            
            # Check if barely changing
            changes = np.diff(recent_values)
            max_change = np.max(np.abs(changes))
            avg_value = np.mean(np.abs(recent_values))
            
            if max_change < threshold * max(avg_value, 1.0):
                stagnant.append(var_name)
        
        return stagnant
    
    def _generate_suggestions(self, diagnostics: Dict) -> List[str]:
        """
        Generate suggestions based on diagnostics
        """
        suggestions = []
        
        # Convergence rate
        rate = diagnostics.get('convergence_rate')
        if rate is not None:
            if rate > 1.0:
                suggestions.append("System is diverging. Try:")
                suggestions.append("  - Reduce relaxation (use smaller steps)")
                suggestions.append("  - Check initial guess (may be far from solution)")
                suggestions.append("  - Verify equations are formulated correctly")
            elif rate > 0.95:
                suggestions.append("Slow convergence. Try:")
                suggestions.append("  - Improve initial guess")
                suggestions.append("  - Use better preconditioning")
                suggestions.append("  - Check for ill-conditioning")
        
        # Problematic equations
        if 'problematic_equations' in diagnostics:
            suggestions.append("\nProblematic equations found:")
            for eq_name, res in diagnostics['problematic_equations'][:3]:
                suggestions.append(f"  - {eq_name}: residual = {res:.6e}")
            suggestions.append("  Check these equations for errors or scaling issues")
        
        # Oscillating variables
        if 'oscillating_variables' in diagnostics:
            suggestions.append("\nOscillating variables detected:")
            for var_name in diagnostics['oscillating_variables'][:3]:
                suggestions.append(f"  - {var_name}")
            suggestions.append("  Try increasing damping/relaxation factor")
        
        # Stagnant variables
        if 'stagnant_variables' in diagnostics:
            suggestions.append("\nStagnant variables detected:")
            for var_name in diagnostics['stagnant_variables'][:3]:
                suggestions.append(f"  - {var_name}")
            suggestions.append("  These may be over-constrained or insensitive")
        
        # Ill-conditioning
        if diagnostics.get('warning_ill_conditioned'):
            suggestions.append("\nPossible ill-conditioning detected:")
            suggestions.append("  - Check variable scaling (normalize variables)")
            suggestions.append("  - Use better numerical methods (e.g., trust region)")
            suggestions.append("  - Add regularization if needed")
        
        return suggestions
    
    def plot_convergence_history(self, max_equations: int = 10, max_variables: int = 10):
        """
        Plot convergence history
        """
        if not self.iteration_data:
            logger.warning("No iteration data to plot")
            return
        
        fig = plt.figure(figsize=(16, 10))
        
        # 1. Overall residual norm
        ax1 = plt.subplot(2, 3, 1)
        iterations = [data['iteration'] for data in self.iteration_data]
        residual_norms = [data['residual_norm'] for data in self.iteration_data]
        ax1.semilogy(iterations, residual_norms, 'o-', linewidth=2)
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('Residual Norm')
        ax1.set_title('Overall Convergence')
        ax1.grid(True, alpha=0.3)
        
        # 2. Individual residuals
        ax2 = plt.subplot(2, 3, 2)
        n_eqs = min(max_equations, len(self.equation_names))
        for i in range(n_eqs):
            eq_name = self.equation_names[i] if i < len(self.equation_names) else f"eq_{i}"
            if eq_name in self.residual_history:
                iters, vals = zip(*self.residual_history[eq_name])
                ax2.semilogy(iters, np.abs(vals), 'o-', label=eq_name, markersize=3)
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('|Residual|')
        ax2.set_title('Individual Equation Residuals')
        ax2.legend(fontsize=8, loc='best')
        ax2.grid(True, alpha=0.3)
        
        # 3. Variable evolution
        ax3 = plt.subplot(2, 3, 3)
        n_vars = min(max_variables, len(self.variable_names))
        for i in range(n_vars):
            var_name = self.variable_names[i] if i < len(self.variable_names) else f"var_{i}"
            if var_name in self.variable_history:
                iters, vals = zip(*self.variable_history[var_name])
                ax3.plot(iters, vals, 'o-', label=var_name, markersize=3)
        ax3.set_xlabel('Iteration')
        ax3.set_ylabel('Variable Value')
        ax3.set_title('Variable Evolution')
        ax3.legend(fontsize=8, loc='best')
        ax3.grid(True, alpha=0.3)
        
        # 4. Residual reduction rate
        ax4 = plt.subplot(2, 3, 4)
        if len(residual_norms) > 1:
            rates = [residual_norms[i]/residual_norms[i-1] if residual_norms[i-1] > 1e-15 else 0
                    for i in range(1, len(residual_norms))]
            ax4.plot(iterations[1:], rates, 'o-')
            ax4.axhline(1.0, color='r', linestyle='--', label='Divergence threshold')
            ax4.axhline(0.9, color='orange', linestyle='--', label='Slow convergence')
            ax4.set_xlabel('Iteration')
            ax4.set_ylabel('Reduction Rate')
            ax4.set_title('Convergence Rate')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        # 5. Largest residuals
        ax5 = plt.subplot(2, 3, 5)
        final_residuals = np.abs(self.iteration_data[-1]['residuals'])
        sorted_indices = np.argsort(final_residuals)[::-1]
        top_n = min(10, len(final_residuals))
        top_indices = sorted_indices[:top_n]
        top_values = final_residuals[top_indices]
        top_names = [self.equation_names[i] if i < len(self.equation_names) else f"eq_{i}" 
                    for i in top_indices]
        ax5.barh(range(top_n), top_values)
        ax5.set_yticks(range(top_n))
        ax5.set_yticklabels(top_names, fontsize=8)
        ax5.set_xlabel('|Residual|')
        ax5.set_title('Largest Final Residuals')
        ax5.set_xscale('log')
        ax5.grid(True, alpha=0.3, axis='x')
        
        # 6. Variable change per iteration
        ax6 = plt.subplot(2, 3, 6)
        if len(self.iteration_data) > 1:
            var_changes = []
            for i in range(1, len(self.iteration_data)):
                change = np.linalg.norm(self.iteration_data[i]['variables'] - 
                                       self.iteration_data[i-1]['variables'])
                var_changes.append(change)
            ax6.semilogy(iterations[1:], var_changes, 'o-')
            ax6.set_xlabel('Iteration')
            ax6.set_ylabel('||Δx||')
            ax6.set_title('Variable Change per Iteration')
            ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def export_report(self, filename: str = 'convergence_report.txt'):
        """
        Export detailed convergence report to file
        """
        diagnostics = self.analyze_convergence()
        
        with open(filename, 'w') as f:
            f.write("="*80 + "\n")
            f.write("CONVERGENCE DIAGNOSTICS REPORT\n")
            f.write("="*80 + "\n\n")
            
            # Overall status
            f.write(f"Converged: {diagnostics.get('converged', False)}\n")
            f.write(f"Final Residual: {diagnostics.get('final_residual', 0):.6e}\n")
            f.write(f"Iterations: {diagnostics.get('n_iterations', 0)}\n")
            
            rate = diagnostics.get('convergence_rate')
            if rate is not None:
                f.write(f"Convergence Rate: {rate:.6f}\n")
            
            f.write("\n" + "-"*80 + "\n\n")
            
            # Warnings
            if 'warning' in diagnostics:
                f.write(f"WARNING: {diagnostics['warning']}\n\n")
            
            # Problematic equations
            if 'problematic_equations' in diagnostics:
                f.write("PROBLEMATIC EQUATIONS:\n")
                for eq_name, res in diagnostics['problematic_equations']:
                    f.write(f"  {eq_name}: {res:.6e}\n")
                f.write("\n")
            
            # Oscillating variables
            if 'oscillating_variables' in diagnostics:
                f.write("OSCILLATING VARIABLES:\n")
                for var_name in diagnostics['oscillating_variables']:
                    f.write(f"  {var_name}\n")
                f.write("\n")
            
            # Suggestions
            if 'suggestions' in diagnostics:
                f.write("SUGGESTIONS:\n")
                for suggestion in diagnostics['suggestions']:
                    f.write(f"{suggestion}\n")
            
            f.write("\n" + "="*80 + "\n")
        
        logger.info(f"Convergence report exported to {filename}")
    
    def clear(self):
        """Clear all recorded data"""
        self.iteration_data = []
        self.variable_history = {}
        self.residual_history = {}


# Convenience function for quick diagnostics
def diagnose_convergence(residual_norms: List[float],
                        tolerance: float = 1e-6) -> str:
    """
    Quick convergence diagnosis from residual history
    
    Args:
        residual_norms: List of residual norms from iterations
        tolerance: Convergence tolerance
        
    Returns:
        Diagnostic message
    """
    if not residual_norms:
        return "No data"
    
    final_residual = residual_norms[-1]
    n_iters = len(residual_norms)
    
    # Check convergence
    if final_residual < tolerance:
        return f"✓ Converged in {n_iters} iterations (residual={final_residual:.6e})"
    
    # Check divergence
    if len(residual_norms) > 2:
        recent = residual_norms[-5:]
        if all(recent[i] > recent[i-1] for i in range(1, len(recent))):
            return f"✗ Diverging! ({n_iters} iterations, residual={final_residual:.6e})"
    
    # Check stagnation
    if len(residual_norms) > 5:
        recent = residual_norms[-5:]
        if max(recent) - min(recent) < 0.01 * np.mean(recent):
            return f"⚠ Stagnated at residual={final_residual:.6e} after {n_iters} iterations"
    
    # Slow convergence
    return f"⚠ Slow convergence: {n_iters} iterations, residual={final_residual:.6e}"


if __name__ == "__main__":
    # Test: Simulate convergence data
    np.random.seed(42)
    
    diag = ConvergenceDiagnostics()
    
    # Simulate iterations with oscillating behavior
    n_vars = 5
    n_eqs = 5
    
    var_names = [f"P_{i}" for i in range(3)] + [f"T_{i}" for i in range(2)]
    eq_names = [f"mass_bal_{i}" for i in range(3)] + [f"energy_bal_{i}" for i in range(2)]
    
    x = np.array([10.0, 8.0, 5.0, 300.0, 295.0])
    
    for iteration in range(20):
        # Simulate convergence with some oscillation
        noise = 0.5 * np.exp(-iteration/5) * np.random.randn(n_vars)
        x += noise
        
        # Compute mock residuals
        residuals = 0.1 * np.exp(-iteration/3) * np.random.randn(n_eqs)
        
        # Record
        diag.record_iteration(iteration, x, residuals, var_names, eq_names)
    
    # Analyze
    diagnostics = diag.analyze_convergence()
    
    print("="*80)
    print("CONVERGENCE DIAGNOSTICS")
    print("="*80)
    for key, value in diagnostics.items():
        if key != 'suggestions':
            print(f"{key}: {value}")
    
    print("\nSUGGESTIONS:")
    for suggestion in diagnostics.get('suggestions', []):
        print(suggestion)
    
    # Plot
    diag.plot_convergence_history()
    
    # Export
    diag.export_report('test_convergence_report.txt')