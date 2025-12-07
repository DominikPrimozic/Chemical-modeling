# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 15:41:54 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# Define box size (a) and circle diameter
  # Adjust as needed
d = 1  # Circle diameter
den = 0.012 * 35
step = 0.04 * 35
N = 100
a = (100/den)**(1/2)
b=a
a+=a/10

# Read the coordinates from the file
data = []
filepath = f"output/random_walk/configuration/{N}_{den:.6f}_0.txt"
coords = np.loadtxt(filepath, unpack=True)

# Function to calculate distance between two points
def distance(p1, p2):
    return np.linalg.norm(p1 - p2)

# Check for overlap
overlap_found = False
for i in range(len(coords)):
    for j in range(i + 1, len(coords)):  # Compare each pair only once
        dist = distance(coords[i], coords[j])
        if dist < d:  # If distance is less than diameter, there is overlap
            print(f"Overlap detected between circles at {coords[i]} and {coords[j]}")
            overlap_found = True

if not overlap_found:
    print("No overlap detected.")

# Create circle mesh
def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# Load the walk data
filepath = f"output/random_walk/path/{N}_{den:.6f}_{step:.6f}_5.txt"
walk = np.loadtxt(filepath, unpack=True)

# Create figure
fig = go.Figure()

# Add circles
for cx, cy in zip(coords[0, :], coords[1, :]):
    # Only plot circles that are within the defined box (i.e., within [-a/2, a/2])
    if abs(cx) <= a/2 and abs(cy) <= a/2:
        x, y = create_circle(cx, cy, r=d/2)
        fig.add_trace(go.Scatter(
            x=x, y=y,
            fill="toself",
            mode="lines",  # Draw the boundary of each circle
            line=dict(color="blue"),
            opacity=0.8,
            showlegend=False  # Remove legend for this trace
        ))

# Add boundary lines at -a/2 and +a/2 in both directions
boundary_lines = [
    # Horizontal lines
    go.Scatter(x=[-b/2, b/2], y=[-b/2, -b/2], mode="lines", line=dict(color="black", dash="dash"), showlegend=False),  # bottom line
    go.Scatter(x=[-b/2, b/2], y=[b/2, b/2], mode="lines", line=dict(color="black", dash="dash"), showlegend=False),  # top line

    # Vertical lines
    go.Scatter(x=[-b/2, -b/2], y=[-b/2, b/2], mode="lines", line=dict(color="black", dash="dash"), showlegend=False),  # left line
    go.Scatter(x=[b/2, b/2], y=[-b/2, b/2], mode="lines", line=dict(color="black", dash="dash"), showlegend=False)   # right line
]

# Add boundary lines to the figure
for line in boundary_lines:
    fig.add_trace(line)

# Create frames for the animation, using separate traces for each frame
frames = []
for i in range(1, len(walk[0])):
    frames.append(go.Frame(
    data=[
        # Redraw all static circles in every frame
        *[
            go.Scatter(
                x=create_circle(cx, cy, r=d/2)[0],
                y=create_circle(cx, cy, r=d/2)[1],
                fill="toself",
                mode="lines",
                line=dict(color="blue"),
                opacity=0.8,
                showlegend=False
            ) for cx, cy in zip(coords[0, :], coords[1, :])
        ],
        # Add the animated red marker
        go.Scatter(
            x=[walk[0, i-1]],
            y=[walk[1, i-1]],
            mode="markers",
            marker=dict(color="red", size=6)
        )
    ],
    name=f"frame_{i}"
))

# Set up the layout
fig.update_layout(
    xaxis=dict(
        range=[-a/2, a/2],  # Set the x-axis range from -a/2 to a/2
        showgrid=True,  # Show grid lines
        zeroline=False,  # Remove zero line (black line in the center)
        tickvals=np.linspace(-a/2, a/2, 5),  # Set ticks at even intervals
        ticktext=[f"{i:.1f}" for i in np.linspace(-a/2, a/2, 5)],  # Customize tick labels
        tickangle=0,  # Angle of the tick labels
        ticks="outside",  # Draw ticks outside the plot area
        tickwidth=1,  # Width of the ticks
        ticklen=8,  # Length of the ticks
        title="x",  # Label for x-axis
    ),
    yaxis=dict(
        range=[-a/2, a/2],  # Set the y-axis range from -a/2 to a/2
        showgrid=True,  # Show grid lines
        zeroline=False,  # Remove zero line (black line in the center)
        tickvals=np.linspace(-a/2, a/2, 5),  # Set ticks at even intervals
        ticktext=[f"{i:.1f}" for i in np.linspace(-a/2, a/2, 5)],  # Customize tick labels
        tickangle=0,  # Angle of the tick labels
        ticks="outside",  # Draw ticks outside the plot area
        tickwidth=1,  # Width of the ticks
        ticklen=8,  # Length of the ticks
        title="y",  # Label for y-axis
    ),
    title=f"{N}_{den:.6f}_{step:.6f}",  # Title of the plot
    width=600,  # Set fixed width for the plot
    height=600,  # Set fixed height for the plot
    plot_bgcolor="white",  # Set background color of the plot
    paper_bgcolor="white",  # Set background color of the paper
    updatemenus=[dict(
        type="buttons",
        showactive=False,
        buttons=[dict(
            label="Play",
            method="animate",
            args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True)]
        )]
    )]
)

# Add the frames to the layout
fig.frames = frames

# Open the animation in a browser
pio.show(fig, renderer='browser')

fig.write_html(f"{N}_{den:.6f}_{step:.6f}.html")
