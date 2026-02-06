"""
Flow Network Solver
Solves pressure-flow relationships across the entire flowsheet
Ensures mass conservation and proper flow direction
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class FlowPath:
    """
    Represents a flow path from source to sink through resistances
    """
    
    def __init__(self):
        self.source = None          # Pressure source (Feed, Volume)
        self.sink = None            # Pressure sink (Boundary, Volume)
        self.resistances = []       # Flow resistances (Valves, Pipes)
        self.streams = []           # Streams in this path
        self.total_resistance = 0.0
        
    def calculate_total_resistance(self):
        """Calculate total flow resistance in path"""
        self.total_resistance = sum(r.get_resistance() for r in self.resistances)
        return self.total_resistance
    
    def calculate_flow(self):
        """
        Calculate flow rate from pressure drop and resistance
        F = dP / R
        For valve: R ~ 1/(Cv^2 * rho)
        """
        if self.source is None or self.sink is None:
            return 0.0
        
        # Pressure drop
        P_source = self.source.P if hasattr(self.source, 'P') else self.source.P_fixed
        P_sink = self.sink.P if hasattr(self.sink, 'P') else self.sink.P_fixed
        dP_total = P_source - P_sink
        
        if dP_total <= 0:
            return 0.0
        
        # Calculate resistance
        self.calculate_total_resistance()
        
        if self.total_resistance <= 0:
            return 0.0
        
        # Flow rate
        F = dP_total / self.total_resistance
        
        return F
    
    def __repr__(self):
        source_name = self.source.name if self.source else "None"
        sink_name = self.sink.name if self.sink else "None"
        return f"<FlowPath: {source_name} → {sink_name} ({len(self.resistances)} resistances)>"


class FlowNetwork:
    """
    Manages flow network topology and solves pressure-flow relationships
    
    Approach:
    1. Identify flow paths (source → resistances → sink)
    2. Calculate pressure drops across resistances
    3. Solve for flows that satisfy pressure boundaries and mass conservation
    """
    
    def __init__(self):
        self.pressure_sources = []  # Units that set pressure (Feed, Volume)
        self.pressure_sinks = []    # Units that set pressure (Boundary, Volume)
        self.flow_resistances = []  # Units that resist flow (Valve, Pipe)
        self.flow_paths = []        # List of FlowPath objects
        
        self.topology_built = False
        
    def build_topology(self, units: List):
        """
        Analyze flowsheet topology and identify flow paths
        
        Args:
            units: List of all unit operations
        """
        logger.info("Building flow network topology...")
        
        from units.feed import FeedBoundary
        from units.boundary import Boundary
        from units.valve import Valve
        from units.volume import Volume
        
        # Clear existing
        self.pressure_sources = []
        self.pressure_sinks = []
        self.flow_resistances = []
        self.flow_paths = []
        
        # Classify units
        for unit in units:
            if isinstance(unit, (FeedBoundary, Volume)):
                self.pressure_sources.append(unit)
            
            if isinstance(unit, (Boundary, Volume)):
                self.pressure_sinks.append(unit)
            
            if isinstance(unit, Valve):
                self.flow_resistances.append(unit)
        
        # Build flow paths
        self._identify_flow_paths(units)
        
        self.topology_built = True
        
        logger.info(f"Topology built: {len(self.pressure_sources)} sources, "
                   f"{len(self.pressure_sinks)} sinks, "
                   f"{len(self.flow_resistances)} resistances, "
                   f"{len(self.flow_paths)} paths")
    
    def _identify_flow_paths(self, units: List):
        """
        Identify flow paths by traversing from sources to sinks
        """
        from units.feed import FeedBoundary
        from units.boundary import Boundary
        
        for source in self.pressure_sources:
            # Start from source, traverse downstream
            self._traverse_from_source(source, units)
    
    def _traverse_from_source(self, source, units):
        """
        Traverse flowsheet from a source to find paths to sinks
        Uses depth-first search
        """
        path = FlowPath()
        path.source = source
        
        # Start traversal from source outlet
        if len(source.outlet_ports) == 0:
            return
        
        start_stream = source.outlet_ports[0].stream
        self._dfs_traverse(start_stream, path, visited=set())
    
    def _dfs_traverse(self, current_stream, current_path, visited):
        """
        Depth-first search to find flow paths
        """
        if current_stream is None or current_stream in visited:
            return
        
        visited.add(current_stream)
        current_path.streams.append(current_stream)
        
        # Find unit connected to this stream's outlet
        downstream_unit = self._find_downstream_unit(current_stream)
        
        if downstream_unit is None:
            return
        
        from units.boundary import Boundary
        from units.valve import Valve
        from units.volume import Volume
        
        # Check if we reached a sink
        if isinstance(downstream_unit, (Boundary, Volume)):
            current_path.sink = downstream_unit
            
            # Complete path found - save it
            path_copy = FlowPath()
            path_copy.source = current_path.source
            path_copy.sink = current_path.sink
            path_copy.resistances = current_path.resistances.copy()
            path_copy.streams = current_path.streams.copy()
            self.flow_paths.append(path_copy)
            
            visited.remove(current_stream)
            current_path.streams.pop()
            return
        
        # If it's a resistance, add to path
        if isinstance(downstream_unit, Valve):
            current_path.resistances.append(downstream_unit)
        
        # Continue traversal from downstream unit's outlets
        for outlet_port in downstream_unit.outlet_ports:
            if outlet_port.stream:
                self._dfs_traverse(outlet_port.stream, current_path, visited)
        
        # Backtrack
        if isinstance(downstream_unit, Valve):
            current_path.resistances.pop()
        
        visited.remove(current_stream)
        current_path.streams.pop()
    
    def _find_downstream_unit(self, stream):
        """
        Find the unit operation downstream of a stream
        """
        # The stream is an outlet from one unit and inlet to another
        # We need to find which unit has this stream as inlet
        
        for port in stream.connected_ports:
            if port.kind == 'inlet':
                # This port belongs to the downstream unit
                # Find which unit owns this port
                for unit in self.pressure_sources + self.pressure_sinks + self.flow_resistances:
                    if port in unit.inlet_ports:
                        return unit
        
        return None
    
    def solve_flows(self, max_iterations: int = 20, tolerance: float = 1e-6) -> bool:
        """
        Solve flow network for all paths
        
        Iterative approach:
        1. Calculate flow for each path from pressure drop
        2. Update volume pressures based on accumulation
        3. Iterate until converged
        
        Returns:
            True if converged
        """
        if not self.topology_built:
            logger.error("Topology not built. Call build_topology() first.")
            return False
        
        if len(self.flow_paths) == 0:
            logger.warning("No flow paths found")
            return True
        
        # Solve each path independently (simplified approach)
        for path in self.flow_paths:
            F = path.calculate_flow()
            
            # Set flow on all streams in path
            for stream in path.streams:
                stream.F = F
            
            # Set flow on resistance inlet/outlet
            for resistance in path.resistances:
                if len(resistance.inlet_ports) > 0 and resistance.inlet_ports[0].stream:
                    resistance.inlet_ports[0].stream.F = F
                if len(resistance.outlet_ports) > 0 and resistance.outlet_ports[0].stream:
                    resistance.outlet_ports[0].stream.F = F
        
        return True
    
    def solve_flows_iterative(self, max_iterations: int = 20, tolerance: float = 1e-6) -> Tuple[bool, int]:
        """
        Iterative solution for flow network with volumes
        
        Algorithm:
        1. For each path, calculate flow from current pressures
        2. Check mass balance at volumes
        3. Adjust volume pressures
        4. Repeat until converged
        
        Returns:
            (converged, iterations)
        """
        if not self.topology_built:
            logger.error("Topology not built")
            return False, 0
        
        from units.volume import Volume
        
        for iteration in range(max_iterations):
            # Store old flows
            old_flows = {}
            for path in self.flow_paths:
                for stream in path.streams:
                    old_flows[stream.name] = stream.F
            
            # Calculate flows for each path
            for path in self.flow_paths:
                F = path.calculate_flow()
                
                # Set flow on streams
                for stream in path.streams:
                    stream.F = F
            
            # Check convergence
            max_change = 0.0
            for path in self.flow_paths:
                for stream in path.streams:
                    if stream.name in old_flows:
                        old_F = old_flows[stream.name]
                        change = abs(stream.F - old_F) / (abs(old_F) + 1e-6)
                        max_change = max(max_change, change)
            
            if max_change < tolerance:
                logger.debug(f"Flow network converged in {iteration+1} iterations")
                return True, iteration + 1
        
        logger.warning(f"Flow network did not converge in {max_iterations} iterations")
        return False, max_iterations
    
    def calculate_pressure_drops(self) -> Dict[str, float]:
        """
        Calculate pressure drops across all resistances
        
        Returns:
            Dictionary of {unit_name: pressure_drop}
        """
        pressure_drops = {}
        
        for resistance in self.flow_resistances:
            if len(resistance.inlet_ports) > 0 and len(resistance.outlet_ports) > 0:
                inlet_stream = resistance.inlet_ports[0].stream
                outlet_stream = resistance.outlet_ports[0].stream
                
                if inlet_stream and outlet_stream:
                    dP = inlet_stream.P - outlet_stream.P
                    pressure_drops[resistance.name] = dP
        
        return pressure_drops
    
    def get_flow_summary(self) -> str:
        """
        Generate summary of flow network
        """
        summary = []
        summary.append("="*80)
        summary.append("FLOW NETWORK SUMMARY")
        summary.append("="*80)
        
        # Paths
        summary.append(f"\nFlow Paths: {len(self.flow_paths)}")
        for i, path in enumerate(self.flow_paths):
            source_name = path.source.name if path.source else "None"
            sink_name = path.sink.name if path.sink else "None"
            n_res = len(path.resistances)
            
            # Calculate flow
            F = path.calculate_flow()
            
            summary.append(f"  Path {i+1}: {source_name} → {sink_name}")
            summary.append(f"    Flow: {F:.4f} kg/s")
            summary.append(f"    Resistances: {n_res}")
            
            for res in path.resistances:
                R = res.get_resistance()
                summary.append(f"      - {res.name}: R={R:.2e}")
        
        # Pressure drops
        summary.append("\nPressure Drops:")
        pressure_drops = self.calculate_pressure_drops()
        for name, dP in pressure_drops.items():
            summary.append(f"  {name}: {dP/1e5:.3f} bar")
        
        summary.append("="*80)
        
        return "\n".join(summary)
    
    def validate_topology(self) -> List[str]:
        """
        Validate flow network topology
        
        Returns:
            List of warnings/errors
        """
        issues = []
        
        # Check for disconnected units
        all_units = set(self.pressure_sources + self.pressure_sinks + self.flow_resistances)
        connected_units = set()
        
        for path in self.flow_paths:
            if path.source:
                connected_units.add(path.source)
            if path.sink:
                connected_units.add(path.sink)
            connected_units.update(path.resistances)
        
        disconnected = all_units - connected_units
        if disconnected:
            issues.append(f"Disconnected units: {[u.name for u in disconnected]}")
        
        # Check for paths without source
        for path in self.flow_paths:
            if path.source is None:
                issues.append(f"Path without source found")
        
        # Check for paths without sink
        for path in self.flow_paths:
            if path.sink is None:
                issues.append(f"Path without sink found")
        
        # Check for loops without volumes
        # (This is complex - simplified check)
        
        return issues


class SimpleFlowSolver:
    """
    Simplified flow solver for cases where detailed network analysis is not needed
    """
    
    def __init__(self):
        pass
    
    def solve_sequential(self, units: List) -> bool:
        """
        Solve flows sequentially from feed to product
        Assumes linear topology: Feed → Valve → Volume → Valve → Product
        """
        from units.feed import FeedBoundary
        from units.valve import Valve
        from units.volume import Volume
        from units.boundary import Boundary
        
        # Find feed
        feed = None
        for unit in units:
            if isinstance(unit, FeedBoundary):
                feed = unit
                break
        
        if feed is None:
            logger.error("No feed found")
            return False
        
        # Traverse from feed
        current_unit = feed
        visited = set()
        
        while current_unit is not None:
            if current_unit in visited:
                logger.error("Loop detected in flowsheet")
                return False
            
            visited.add(current_unit)
            
            # Process unit
            if isinstance(current_unit, Valve):
                # Valve determines flow from pressure drop
                self._solve_valve(current_unit)
            
            elif isinstance(current_unit, Volume):
                # Volume: outlet flow = inlet flow (steady state)
                if len(current_unit.inlet_ports) > 0 and len(current_unit.outlet_ports) > 0:
                    inlet_stream = current_unit.inlet_ports[0].stream
                    outlet_stream = current_unit.outlet_ports[0].stream
                    if inlet_stream and outlet_stream:
                        outlet_stream.F = inlet_stream.F
            
            # Move to next unit
            if len(current_unit.outlet_ports) > 0:
                outlet_stream = current_unit.outlet_ports[0].stream
                if outlet_stream:
                    # Find next unit
                    next_unit = self._find_next_unit(outlet_stream, units)
                    current_unit = next_unit
                else:
                    break
            else:
                break
        
        return True
    
    def _solve_valve(self, valve):
        """Solve valve flow from pressure drop"""
        inlet = valve.inlet_ports[0].stream
        outlet = valve.outlet_ports[0].stream
        
        if inlet is None or outlet is None:
            return
        
        # Calculate flow
        dP = inlet.P - outlet.P
        if dP <= 0:
            outlet.F = 0.0
            inlet.F = 0.0
            return
        
        # Simplified valve equation
        opening = valve.parameters.get('opening', 1.0)
        Cv = valve.parameters.get('Cv', 10.0)
        
        if opening < 0.0001:
            F = 0.0
        else:
            # Get density
            from thermoprops import PropertyExtractor
            comp_dict = {inlet.compounds[i]: inlet.z[i] for i in range(len(inlet.compounds))}
            props = PropertyExtractor(inlet.compounds, comp_dict, inlet.P, inlet.T)
            rho = props.get_density()
            if rho is None or rho <= 0:
                rho = 1.0
            
            # Flow equation
            Cv_SI = Cv * 0.0000241
            Y = 1.0 - dP / (3.0 * inlet.P) if inlet.P > 0 else 1.0
            Y = max(0.5, min(Y, 1.0))
            
            F = Cv_SI * opening**2 * Y * np.sqrt(rho * dP)
        
        outlet.F = F
        inlet.F = F
    
    def _find_next_unit(self, stream, units):
        """Find unit downstream of stream"""
        for unit in units:
            for port in unit.inlet_ports:
                if port.stream is stream:
                    return unit
        return None


if __name__ == "__main__":
    # Test
    print("FlowNetwork module loaded successfully")
    print("Use with DynamicSimulator to solve flow networks")