def ask_and_plot(simulator, time_points):
    import tkinter as tk
    from tkinter import ttk
    import matplotlib.pyplot as plt
    import numpy as np

    # Gather all objects: streams and units
    all_objs = {}
    # Streams
    for name, stream in simulator.streams.items():
        vars_list = [attr for attr in dir(stream) if attr.endswith('_history')]
        all_objs[f"Stream: {name}"] = (stream, vars_list)
    # Units
    for unit in simulator.units:
        vars_list = [attr for attr in dir(unit) if attr.endswith('_history')]
        all_objs[f"Unit: {unit.name}"] = (unit, vars_list)

    obj_names = list(all_objs.keys())

    root = tk.Tk()
    root.title("Plot Simulation Results")

    tk.Label(root, text="Select Object:").grid(row=0, column=0, sticky="w")
    obj_var = tk.StringVar(value=obj_names[0])
    obj_menu = ttk.Combobox(root, textvariable=obj_var, values=obj_names, state="readonly")
    obj_menu.grid(row=0, column=1)

    tk.Label(root, text="Select Variable:").grid(row=1, column=0, sticky="w")
    var_var = tk.StringVar(value=all_objs[obj_names[0]][1][0])
    var_menu = ttk.Combobox(root, textvariable=var_var, values=all_objs[obj_names[0]][1], state="readonly")
    var_menu.grid(row=1, column=1)

    def update_var_menu(*args):
        vars_for_obj = all_objs[obj_var.get()][1]
        var_menu['values'] = vars_for_obj
        var_var.set(vars_for_obj[0])
    obj_var.trace('w', update_var_menu)

    def plot_selected():
        obj, _ = all_objs[obj_var.get()]
        y = getattr(obj, var_var.get())
        times = np.array(time_points)
        min_len = min(len(times), len(y))
        times = times[:min_len]
        y = y[:min_len]
        plt.figure(figsize=(8, 5))
        plt.plot(times, y, '-o')
        plt.xlabel("Time (s)")
        plt.ylabel(var_var.get())
        plt.title(f"{obj_var.get()} - {var_var.get()} vs Time")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    plot_btn = tk.Button(root, text="Plot", command=plot_selected)
    plot_btn.grid(row=2, column=0, columnspan=2, pady=10)

    root.mainloop()