# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 18:02:02 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import re
import os

# Function to read simulation parameters (unchanged)
def read_simulation_params(filepath):
    with open(filepath, 'r') as f:
        N = int(f.readline().split('//')[0].strip())
        density = float(f.readline().split('//')[0].strip())
        f.readline()
        D = int(f.readline().split('//')[0].strip())
    return N, density, D

# Read parameters
param_file = "input/MDLJ/input.txt"
N, density, D = read_simulation_params(param_file)
a = (N / density) ** (1 / D)
d = 1  # Sphere diameter

# Function to read particle data (unchanged)
def read_particles_binary(filename, N, D):
    data = np.fromfile(filename, dtype=np.float64)
    assert data.size == N * D, f"Expected {N * D} elements, got {data.size}"
    return data.reshape((N, D))

# Automatically gather all .bin files in the directory
directory = "output/MDLJ/T5/pos/"  # Your directory containing .bin files
filenames = sorted(
    [f for f in os.listdir(directory) if f.endswith('.bin')],
    key=lambda x: float(re.search(r'(\d+\.\d+)\.bin$', x).group(1))  # Sort numerically
)

# Extract the step part from the filenames
steps = [re.search(r'(\d+\.\d+)\.bin$', filename).group(1) for filename in filenames]

# Create circle mesh (2D) for plotting
def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# Create figure
fig = go.Figure()

# Function to add circles to the figure
def add_circles_for_configuration(coords):
    traces = []
    for cx, cy in coords:
        if abs(cx) <= a / 2 and abs(cy) <= a / 2:
            x, y = create_circle(cx, cy, r=d / 2)
            traces.append(go.Scatter(
                x=x, y=y,
                fill="none",  # No fill, just outline
                mode="lines",
                line=dict(color='blue', dash='dot'),  # Dashed outline
                opacity=0.6,
                showlegend=False
            ))

            # Add a smaller dot at the center of each circle
            traces.append(go.Scatter(
                x=[cx], y=[cy],
                mode="markers",
                marker=dict(size=2, color="red", symbol="circle"),
                showlegend=False
            ))

    return traces

# Add traces for the bounding box (edges)
def add_bounding_box():
    box_traces = [
        # Left vertical line (x = -a/2)
        go.Scatter(
            x=[-a / 2, -a / 2],
            y=[-a / 2, a / 2],
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False
        ),
        # Right vertical line (x = a/2)
        go.Scatter(
            x=[a / 2, a / 2],
            y=[-a / 2, a / 2],
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False
        ),
        # Bottom horizontal line (y = -a/2)
        go.Scatter(
            x=[-a / 2, a / 2],
            y=[-a / 2, -a / 2],
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False
        ),
        # Top horizontal line (y = a/2)
        go.Scatter(
            x=[-a / 2, a / 2],
            y=[a / 2, a / 2],
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False
        )
    ]
    return box_traces

# Function to update the layout of the figure with the slider and fixed aspect ratio
def update_layout_with_slider():
    fig.update_layout(
        sliders=[{
            'yanchor': 'top', 'xanchor': 'center',  # Corrected xanchor to center
            'currentvalue': {
                'visible': True,
                'prefix': 'Snapshot at ',
                'font': {'size': 20},
                'offset': 10  # Optional: Adds space between the value and the slider
            },
            'steps': [
                {
                    'args': [
                        [step],  # This will be used to update the current snapshot
                        {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate'}
                    ],
                    'label': step,
                    'method': 'animate'
                } for step in steps
            ]
        }],
        # Add aspect ratio constraint to avoid stretching
        xaxis=dict(scaleanchor="y", showgrid=True, zeroline=True),  # This ensures that the x-axis is scaled equally with the y-axis
        yaxis=dict(scaleanchor="x", showgrid=True, zeroline=True),  # This ensures that the y-axis is scaled equally with the x-axis
        plot_bgcolor='white',  # White background for the plot area
        paper_bgcolor='lightgray',  # Light gray background for the paper area
        xaxis_title='X-axis',  # Title for X-axis
        yaxis_title='Y-axis',  # Title for Y-axis
        margin=dict(t=40, b=40, l=40, r=40),  # Adjust margins around the plot
    )

# Add the initial configuration to the plot
coords = read_particles_binary(os.path.join(directory, filenames[0]), N, D)
traces = add_circles_for_configuration(coords)
for trace in traces:
    fig.add_trace(trace)

# Add the bounding box lines
box_traces = add_bounding_box()
for box_trace in box_traces:
    fig.add_trace(box_trace)

# Update layout and sliders
update_layout_with_slider()

# Add the frames (each file represents a snapshot at a given time step)
frames = []
for idx, filename in enumerate(filenames):
    coords = read_particles_binary(os.path.join(directory, filename), N, D)
    frame_traces = add_circles_for_configuration(coords)
    frames.append(go.Frame(data=frame_traces, name=steps[idx]))

fig.frames = frames

# Show the figure
pio.show(fig, renderer='browser')
fig.write_html("T125.html", auto_open=True)