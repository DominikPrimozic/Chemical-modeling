# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 17:45:54 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import re

def read_simulation_params(filepath):
    with open(filepath, 'r') as f:
        N = int(f.readline().split('//')[0].strip())
        density = float(f.readline().split('//')[0].strip())
        f.readline()
        D = int(f.readline().split('//')[0].strip())
    return N, density, D

# Example usage
param_file = "input/MDLJ/input.txt"
N, density, D = read_simulation_params(param_file)
#density=0.1

a = (N / density) ** (1 / D)
d = 1  # Sphere diameter

def read_particles_binary(filename, N, D):
    data = np.fromfile(filename, dtype=np.float64)
    assert data.size == N * D, f"Expected {N * D} elements, got {data.size}"
    return data.reshape((N, D))

# Example usage
filename = "output/MDLJ/T1/pos/0.000000.bin"

coords = read_particles_binary(filename, N, D)


# Use regular expression to extract 'step_2000'
match = re.search(r'(\d+\.\d+)\.bin$', filename)

if match:
    extracted_string = match.group(1)


# Create circle mesh (2D)
def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# Create figure
fig = go.Figure()

# Add circles with periodic boundary conditions
for cx, cy in coords:
    # Only plot circles that are within the defined box (i.e., within [-a/2, a/2])
    if abs(cx) <= a / 2 and abs(cy) <= a / 2:
        x, y = create_circle(cx, cy, r=d / 2)
        fig.add_trace(go.Scatter(
            x=x, y=y,
            fill="none",  # No fill, just outline
            mode="lines",
            line=dict(color='blue', dash='dot'),  # Dashed outline
            opacity=0.6,
            showlegend=False  # Remove legend for this trace
        ))

        # Add a smaller dot at the center of each circle
        fig.add_trace(go.Scatter(
            x=[cx], y=[cy],  # The center of the circle
            mode="markers",
            marker=dict(size=2, color="red", symbol="circle"),  # Smaller size for the dot
            showlegend=False  # Remove legend for this trace
        ))

        # Replicate the particles into the 4 boxes by shifting their coordinates
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx != 0 or dy != 0:  # Skip the original box
                    new_cx = cx + dx * a
                    new_cy = cy + dy * a
                    x, y = create_circle(new_cx, new_cy, r=d / 2)
                    fig.add_trace(go.Scatter(
                        x=x, y=y,
                        fill="none",  # No fill, just outline
                        mode="lines",
                        line=dict(color='blue', dash='dot'),  # Dashed outline
                        opacity=0.6,
                        showlegend=False  # Remove legend for this trace
                    ))

                    # Add a smaller dot at the center of the replicated circles
                    fig.add_trace(go.Scatter(
                        x=[new_cx], y=[new_cy],  # The center of the replicated circle
                        mode="markers",
                        marker=dict(size=2, color="red", symbol="circle"),  # Smaller size for the dot
                        showlegend=False  # Remove legend for this trace
                    ))

# Add rectangle outlines for each box
box_size = a  # Size of each box
for dx in [-1, 0, 1]:
    for dy in [-1, 0, 1]:
        if dx != 0 or dy != 0:  # Skip the original box
            fig.add_trace(go.Scatter(
                x=[-box_size / 2 + dx * box_size, box_size / 2 + dx * box_size, box_size / 2 + dx * box_size, -box_size / 2 + dx * box_size, -box_size / 2 + dx * box_size],
                y=[-box_size / 2 + dy * box_size, -box_size / 2 + dy * box_size, box_size / 2 + dy * box_size, box_size / 2 + dy * box_size, -box_size / 2 + dy * box_size],
                mode='lines',
                line=dict(color='black', width=2),  # Black solid line for the box outline
                showlegend=False
            ))

# Update layout for a bigger box and add axis labels
fig.update_layout(
    xaxis=dict(
        range=[-a * 1.2, a * 1.2],  # Increase the range to make the box bigger
        fixedrange=True,
        showgrid=True,
        zeroline=False,
        tickvals=np.linspace(-a, a, 5),
        ticktext=[f"{i:.1f}" for i in np.linspace(-a, a, 5)],
        tickangle=0,
        ticks="outside",
        tickwidth=1,
        ticklen=8,
        title="x",
    ),
    yaxis=dict(
        range=[-a * 1.2, a * 1.2],  # Increase the range to make the box bigger
        fixedrange=True,
        showgrid=True,
        zeroline=False,
        tickvals=np.linspace(-a, a, 5),
        ticktext=[f"{i:.1f}" for i in np.linspace(-a, a, 5)],
        tickangle=0,
        ticks="outside",
        tickwidth=1,
        ticklen=8,
        title="y",
    ),
    title=f"Snapshot at {extracted_string}",
    width=700,  # Increase width
    height=700,  # Increase height
    plot_bgcolor="white",
    paper_bgcolor="white",
)

# Open in browser
pio.show(fig, renderer='browser')
