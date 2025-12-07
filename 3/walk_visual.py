# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 21:45:54 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 20:40:24 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# Define box size (a) and circle diameter
a = (100/0.3)**(1/2)  # Adjust as needed
d = 1  # Circle diameter
den=0.012*48
step=0.04*49
N=100
# Read the coordinates from the file
data = []
filepath=f"output/random_walk/configuration/{N}_{den:.6f}_0.txt"
coords=np.loadtxt(filepath,unpack=True)


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

filepath=f"output/random_walk/path/{N}_{den:.6f}_{step:.6f}_5.txt"
walk=np.loadtxt(filepath,unpack=True)



# Create figure
fig = go.Figure()

# Add circles
for cx, cy in zip(coords[0,:],coords[1,:]):
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
      #  fixedrange=True,  # Prevent zooming or panning
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
       # fixedrange=True,  # Prevent zooming or panning
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
    title="Diski v 2D",  # Title of the plot
    width=600,  # Set fixed width for the plot
    height=600,  # Set fixed height for the plot
    plot_bgcolor="white",  # Set background color of the plot
    paper_bgcolor="white",  # Set background color of the paper
)

fig.add_trace(go.Scatter(
    x=walk[0, :], 
    y=walk[1, :],  
    mode="lines+markers",  # Add lines between markers
    marker=dict(color="red", size=4),  
    line=dict(color="blue"),  # Line color
    name="Random Walk"
))
# Open in browser
pio.show(fig, renderer='browser')
