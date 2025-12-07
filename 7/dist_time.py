# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 20:11:40 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import os
import webbrowser

def read_particle_positions(filename, num_particles, dimensions=2):
    """
    Read a single configuration from a binary file.
    Returns: (N x dim) array
    """
    positions = np.zeros((num_particles, dimensions))
    with open(filename, "rb") as f:
        for i in range(num_particles):
            positions[i] = np.fromfile(f, dtype=np.float64, count=dimensions)
    return positions

def compute_gr(positions_list, density, dr=0.01, r_max=None, dim=2):
    """
    Compute radial distribution function g(r).

    Args:
        positions_list: list of (N x dim) arrays
        density: number density (N / volume)
        dr: bin width
        r_max: max radius to compute (default: half box)
        dim: 1, 2, or 3

    Returns:
        r: array of r values
        g_r: g(r) values
    """
    num_configs = len(positions_list)
    N = positions_list[0].shape[0]
    volume = N / density
    box_length = volume ** (1 / dim)

    if r_max is None:
        r_max = box_length / 2.0

    N_bins = int(r_max / dr)
    hist = np.zeros(N_bins)

    for pos in positions_list:
        for i in range(N):
            for j in range(i + 1, N):
                rij = pos[i] - pos[j]
                rij -= box_length * np.round(rij / box_length)  # minimum image
                r = np.linalg.norm(rij)
                bin_index = int(r / dr)
                if bin_index < N_bins:
                    hist[bin_index] += 2

    r_vals = (np.arange(N_bins) + 0.5) * dr
    if dim == 3:
        shell_volume = 4 * np.pi * r_vals**2 * dr
    elif dim == 2:
        shell_volume = 2 * np.pi * r_vals * dr
    elif dim == 1:
        shell_volume = 2 * dr
    else:
        raise ValueError("Unsupported dimension")

    g_r = hist / (num_configs * N * density * shell_volume)
    return r_vals, g_r

def read_simulation_params(filepath):
    with open(filepath, 'r') as f:
        N = int(f.readline().split('//')[0].strip())
        density = float(f.readline().split('//')[0].strip())
        D = int(f.readline().split('//')[0].strip())
    return N, density, D

# -------------------------
# 🔧 USER INPUT
config_folder = "output/MDLJ/2/pos/"  # Your folder with .bin files 
#param_file = "input/LJ.txt"
#N, rho, dim = read_simulation_params(param_file)
N=100
rho=0.4
dim=2 
dr = 0.1                    # Bin width

# -------------------------
# 📂 Read all configs from folder
file_paths = sorted([os.path.join(config_folder, f) for f in os.listdir(config_folder) if f.endswith('.bin')])
positions_list = [read_particle_positions(fp, num_particles=N, dimensions=dim) for fp in file_paths]

# -------------------------
# 🚀 Compute g(r) for each configuration
r_values, gr_values_over_time = [], []

for pos in positions_list:
    r, g_r = compute_gr([pos], density=rho, dr=dr, dim=dim)
    r_values = r  # r doesn't change over time
    gr_values_over_time.append(g_r)

# Convert gr_values_over_time to a numpy array for easier manipulation
gr_values_over_time = np.array(gr_values_over_time)

# -------------------------
# 📈 3D Plot using Plotly
fig = go.Figure(data=[
    go.Surface(
        z=gr_values_over_time,  # g(r) values
        x=np.tile(r_values, (len(gr_values_over_time), 1)),  # r values (same for all time steps)
        y=np.arange(len(gr_values_over_time)),  # Time steps
        colorscale='Viridis',
        colorbar=dict(title='g(r)')
    )
])

# Update layout for the plot
fig.update_layout(
    scene=dict(
        xaxis_title='r',
        yaxis_title='Time Step',
        zaxis_title='g(r)'
    ),
    title='Radial Distribution Function over Time',
    autosize=True
)

# Save the figure as an HTML file
output_filename = "radial_distribution_function.html"
fig.write_html(output_filename)

# Automatically open the HTML file in the default web browser
webbrowser.open(output_filename)
