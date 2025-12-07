import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.offline as pyo
import subprocess

def read_input_file(filename):
    data = {}
    variable_names = [
        "c0", "c_0", "c_h", "D", "v_max", "K_m", "h", "T", "Ny", "Nt"
    ]

    with open(filename, "r") as file:
        index = 0
        for line in file:
            line = line.split("//")[0].strip()  # Remove comments and whitespace
            if line:
                try:
                    data[variable_names[index]] = float(line)
                    index += 1
                except ValueError:
                    print(f"Warning: Could not convert '{line}' to a number.")
                if index >= len(variable_names):
                    break

    return data

input_file = "input/input.txt"
values = read_input_file(input_file)

dy = 1 / values["Ny"]
dt = 1 / values["Nt"]
data_path = f"output/reactor/{values['D']:.6f}_{values['v_max']:.6f}_{values['K_m']:.6f}_{dy:.6f}_{dt:.6f}.txt"

subprocess.run(["run_reaction.exe"], check=True)

data = np.loadtxt(data_path, delimiter="\t")

time = data[:, 0]
concentration = data[:, 1:]

x = np.linspace(0, 1, concentration.shape[1])

# Create figure
fig = go.Figure()

# Add frames for animation
frames = [go.Frame(
    data=[go.Scatter(x=x, y=concentration[i, :], mode='lines', name=f'Time {time[i]:.2f}')],
    name=str(i)
) for i in range(len(time))]

fig.frames = frames

# Add initial trace
fig.add_trace(go.Scatter(x=x, y=concentration[0, :], mode='lines', name='Initial'))

# Set layout with animation settings
fig.update_layout(
    title="Concentration Profile Over Time",
    xaxis_title="y",
    yaxis_title="z",
    yaxis=dict(range=[0, 1]),  # Fixed concentration axis from 0 to 1
    updatemenus=[
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ],
    sliders=[
        {
            "active": 0,
            "yanchor": "top",
            "xanchor": "left",
            "currentvalue": {"font": {"size": 20}, "prefix": "Time: ", "visible": True, "xanchor": "right"},
            "pad": {"b": 10, "t": 50},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [[f.name], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}],
                    "label": f"{time[k]:.2f}",
                    "method": "animate"
                } for k, f in enumerate(fig.frames)
            ]
        }
    ]
)

# Save and show animation
pyo.plot(fig, filename=f"animation_concentration_{values['D']:.6f}_{values['v_max']:.6f}_{values['K_m']:.6f}_{dy:.6f}_{dt:.6f}.html")