# units/unit_operation.py
import os
import re 
import config.sim_settings as sim_settings  # At the top of the file
import builtins
class Port:
    """
    Represents a port on a unit operation.
    Can be 'inlet' or 'outlet', and is connected to a Stream.
    """
    def __init__(self, name, kind, stream=None):
        assert kind in ('inlet', 'outlet')
        self.name = name
        self.kind = kind
        self.stream = stream


    def connect(self, stream):
        self.stream = stream
        stream.connect_port(self)

    def __repr__(self):
        return f"<Port {self.name} ({self.kind})>"

class UnitOperation:
    """
    Base class for all unit operations.
    Supports any number of inlet and outlet ports.
    """
    def __init__(self, name, n_inlets=1, n_outlets=1):
        self.name = name
        self.inlet_ports = [Port(f"{name}_in{i+1}", 'inlet') for i in range(n_inlets)]
        self.outlet_ports = [Port(f"{name}_out{i+1}", 'outlet') for i in range(n_outlets)]
        self.parameters = {}
        # User code slots (None or string)
        self.user_code_before = None
        self.user_code_after = None
        self.arrayed_vars = {}  # For user-defined arrayed variables

    def store_history(self):
        # Store user-defined arrayed variables
        for var, value in self.arrayed_vars.items():
            hist_name = f"{var}_history"
            if not hasattr(self, hist_name):
                setattr(self, hist_name, [])
            getattr(self, hist_name).append(value)

    def load_user_code_from_master_file(self, code_file='user_code.py'):
        """
        Loads user code for this unit from a single master file.
        Looks for blocks marked as:
        # --- <UnitName>_before ---
        # --- <UnitName>_after ---
        """
        if not os.path.isfile(code_file):
            return

        import chardet

        with open(code_file, 'rb') as f:
            raw_data = f.read()

        detected = chardet.detect(raw_data)
        encoding = detected['encoding']

        with open(code_file, 'r', encoding=encoding) as f:
            code_text = f.read()



        # Regex to extract code blocks
        for when in ['before', 'after']:
            pattern = rf"# --- {re.escape(self.name)}_{when} ---\n(.*?)(?=\n# ---|\Z)"
            match = re.search(pattern, code_text, re.DOTALL)
            if match:
                code = match.group(1).strip()
                if when == 'before':
                    self.user_code_before = code if code else None
                else:
                    self.user_code_after = code if code else None
                print(f"Loaded user code for {self.name} {when}: {code[:30]}...")
                print(f"Loaded user code for {self.name} before: {self.user_code_before}")
                print(f"Loaded user code for {self.name} after: {self.user_code_after}")


    def set_user_code(self, code_str, when='before'):
        """Set user code to be executed before or after calculate()."""
        if when == 'before':
            self.user_code_before = code_str
        elif when == 'after':
            self.user_code_after = code_str
        else:
            raise ValueError("when must be 'before' or 'after'")
    def run_user_code(self, when, extra_vars=None):
        code = self.user_code_before if when == 'before' else self.user_code_after
        if code:
            CT = sim_settings.CT
            RT = sim_settings.RT
            # Custom print function
            def print_at_result_interval(*args, **kwargs):
                # Only print if CT is a multiple of RT (within tolerance)
                if abs(CT % RT) < 1e-8 or abs(RT - (CT % RT)) < 1e-8:
                    builtins.print(*args, **kwargs)
            local_vars = {
                'self': self,
                'CT': CT,
                'DT': sim_settings.DT,
                'RT': RT,
                'TOTAL_TIME': sim_settings.TOTAL_TIME,
                'WHEN': when,
                'print': print_at_result_interval  # Override print
            }
            if extra_vars:
                local_vars.update(extra_vars)
            try:
                exec(code, {}, local_vars)
            except Exception as e:
                builtins.print(f"Error in user code for {self.name} {when}: {e}")

    def connect_inlet(self, idx, stream):
        self.inlet_ports[idx].connect(stream)

    def connect_outlet(self, idx, stream):
        self.outlet_ports[idx].connect(stream)

    def get_inlet_streams(self):
        return [port.stream for port in self.inlet_ports]

    def get_outlet_streams(self):
        return [port.stream for port in self.outlet_ports]

    def calculate(self, *args, **kwargs):
        """
        To be implemented by subclasses.
        Should update outlet streams based on inlet streams and unit behavior.
        """
        raise NotImplementedError("Each unit operation must implement its own calculate() method.")

    def __repr__(self):
        return (f"<UnitOperation {self.name}: "
                f"{len(self.inlet_ports)} inlets, {len(self.outlet_ports)} outlets>")
