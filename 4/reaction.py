# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 22:05:43 2025

@author: domin
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.offline as pyo
import subprocess

def read_input_file(filename):
    data = {}
    variable_names = [
        "c0", "c_0", "c_h", "D", "v_max", "K_m", "h", "T", "Ny","Nt"
    ]

    with open(filename, "r") as file:
        index = 0
        for line in file:
            line = line.split("//")[0].strip()  # Remove comments and whitespace
            if line:  # Only process non-empty lines
                try:
                    data[variable_names[index]] = float(line)
                    index += 1
                except ValueError:
                    print(f"Warning: Could not convert '{line}' to a number.")
                if index >= len(variable_names):  # Stop if we have all variables
                    break

    return data

input_file = "input/input.txt"
values = read_input_file(input_file)

dy=1/values["Ny"]
dt=1/values["Nt"]
data_path=f"output/reactor/{values['D']:.6f}_{values['v_max']:.6f}_{values['K_m']:.6f}_{dy:.6f}_{dt:.6f}.txt"

#df = pd.read_csv(data_path, delimiter="\t", header=None)
subprocess.run(["run_reaction.exe"], check=True)

data = np.loadtxt(data_path, delimiter="\t")

time= data[:, 0]
concentration = data[:, 1:]

x=np.linspace(0,1,concentration.shape[1])

X, T = np.meshgrid(x, time)

# Create a 3D surface plot
fig = go.Figure(data=[go.Surface(z=concentration, x=X, y=T)])

# Set labels
fig.update_layout(
    title="Concentration Profile",
    scene=dict(
        xaxis_title="y",
        yaxis_title="t",
        zaxis_title="z"
    )
)

# Render in browser
pyo.plot(fig, filename=f"concentration_{values['D']:.6f}_{values['v_max']:.6f}_{values['K_m']:.6f}_{dy:.6f}_{dt:.6f}.html")
#fig.write_image("3d_concentration.png", scale=2) 

"""
plt.plot(x,concentration[0, :], label="Concentration")

plt.xlabel("y")
plt.ylabel("Concentration")
plt.legend()
plt.show()
"""

"""
time of one conc
plt.plot(time, concentration[:, 4], label="Concentration 5")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
"""