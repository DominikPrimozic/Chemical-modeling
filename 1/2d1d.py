# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 20:40:24 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# Define box size (a) and circle diameter
a = (10/0.5)  # Adjust as needed
d = 1  # Circle diameter

# Read the coordinates from the file
data = []
with open("output/placements/1.10.0.500000.txt", "r") as file:
    for line in file:
        if line.strip() and not line.lower().startswith("packing"):
            values = line.split()
            if len(values) == 1:  # Only process lines with two coordinates (2D)
                data.append([float(values[0]),0])

# Convert list to numpy array
coords = np.array(data)

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

# Create figure
fig = go.Figure()

# Add circles
for cx, cy in coords:
    # Only plot circles that are within the defined box (i.e., within [-a/2, a/2])
    if abs(cx) <= a/2 and abs(cy) <= a/2:
        x, y = create_circle(cx, cy, r=d/2)
        fig.add_trace(go.Scatter(
            x=x, y=y,
            fill="toself",
            mode="lines",
            line=dict(color="blue"),
            opacity=0.8,
            showlegend=False  # Remove legend for this trace
        ))


fig.update_layout(
    xaxis=dict(
        range=[-a/2, a/2],  # Set the x-axis range from -a/2 to a/2
        fixedrange=True,  # Prevent zooming or panning
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
        range=[-0.51, 0.51],  # Set the y-axis range from -a/2 to a/2
        fixedrange=True,  # Prevent zooming or panning
        showgrid=True,  # Show grid lines
        zeroline=False,  # Remove zero line (black line in the center)
        showticklabels=False,  # Remove tick labels
        ticks="",  # Remove ticks
        title="",  # Label for y-axis
    ),
    title="Diski v 1D",  # Title of the plot
    width=600,  # Set fixed width for the plot
    height=600,  # Set fixed height for the plot
    plot_bgcolor="white",  # Set background color of the plot
    paper_bgcolor="white",  # Set background color of the paper
    autosize=False,
    xaxis_scaleanchor="y",  # Lock x-axis scale to y-axis to keep aspect ratio
)


# Open in browser
pio.show(fig, renderer='browser')
