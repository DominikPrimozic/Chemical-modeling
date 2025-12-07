# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 14:20:22 2025

@author: domin
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from moviepy.editor import ImageSequenceClip
import os
# Define box size (a) and circle diameter
d = 1  # Circle diameter
den = 0.012 * 48
step = 0.04 * 49
N = 100
a = (100/den)**(1/2)
b=a
a+=a/10

# Read the coordinates from the file
coords = np.loadtxt(f"output/random_walk/configuration/{N}_{den:.6f}_0.txt", unpack=True)

# Create circle mesh
def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# Load the walk data
walk = np.loadtxt(f"output/random_walk/path/{N}_{den:.6f}_{step:.6f}_5.txt", unpack=True)

# Create figure
fig = go.Figure()

# Add circles
for cx, cy in zip(coords[0, :], coords[1, :]):
    if abs(cx) <= a / 2 and abs(cy) <= a / 2:
        x, y = create_circle(cx, cy, r=d / 2)
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

for line in boundary_lines:
    fig.add_trace(line)

# Generate and save frames for the animation
frame_images = []
for i in range(1, len(walk[0])):
    # Create scatter trace for the current step
    fig.data[0].x = walk[0, :i]  # Update the x data
    fig.data[0].y = walk[1, :i]  # Update the y data

    # Save current frame as an image
    image_filename = f"frame_{i}.png"
    fig.write_image(image_filename)
    frame_images.append(image_filename)

clip = ImageSequenceClip(frame_images, fps=30)
clip.write_videofile("random_walk_animation.mp4", codec="libx264")

# Clean up the generated frame images
for image in frame_images:
    os.remove(image)
