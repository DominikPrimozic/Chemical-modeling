# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 21:05:08 2025

@author: domin
"""


import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go


def read_file(path):
    df = pd.read_csv(path, sep=r'\s+', names=["L", "x", "y"], skiprows=0)
    letters = df["L"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    
    return letters, x, y

def map_labels_to_colors(labels):
    label_colors = {
        'A': 'blue',  # Color for A
        'B': 'green', # Color for B
        'C': 'red'    # Color for C
    }
    return [label_colors[label] for label in labels]

folder_path = 'output/times'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))



fig = go.Figure()
# Read the first file for initialization (frame 0)
letters, x, y = read_file(os.path.join(folder_path, files[0]))  # Include the first frame
color_values = map_labels_to_colors(letters)

# Add initial trace (with just points, no lines)
fig.add_trace(go.Scatter(
    x=x, y=y, mode='markers',  # Use 'markers' for just points
    marker=dict(size=8, color=color_values),  # Color based on labels
    name='Particles'
))

# Create frames for animation (just points)
frames = []
for i, file_name in enumerate(files):
    letters, x, y = read_file(os.path.join(folder_path, file_name))
    color_values = map_labels_to_colors(letters)
    
    # Create frame with just points
    frames.append(go.Frame(
        data=[go.Scatter(
            x=x, y=y, mode='markers',  # Points only (no lines)
            marker=dict(size=8, color=color_values),  # Size and color as before
        )],
        name=str(i),
        layout=go.Layout(
            annotations=[dict(
                x=0.02, y=0.98,
                xref="paper", yref="paper",
                text=f"Frame: {i}",
                showarrow=False,
                font=dict(size=14, color="black"),
                align="left"
            )]
        )
    ))

fig.frames = frames

# Add play button and frame-navigation slider
fig.update_layout(
    title="Animated Reaction with Labels A, B, C",
    xaxis=dict(title="X"),
    yaxis=dict(title="Y"),
    showlegend=True,
    updatemenus=[dict(
        type="buttons",
        showactive=True,
        buttons=[
            dict(
                label="Play",
                method="animate",
                args=[None, dict(frame=dict(duration=200, redraw=True), fromcurrent=True)]
            ),
            dict(
                label="Pause",
                method="animate",
                args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", fromcurrent=True)]
            )
        ]
    )],
    sliders=[{
        "active": 0,  # Set the slider to start at frame 0
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {"prefix": "Frame: ", "font": {"size": 20}, "visible": True},
        "steps": [
            {
                "args": [[str(i)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                "label": str(i),
                "method": "animate"
            }
            for i in range(len(files))  # Include all frames from 0 onward
        ]
    }]
)

fig.show()

fig.write_html("reaction_animation_with_first_frame_visible3.html")