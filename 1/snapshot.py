# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 21:24:49 2025

@author: domin
"""

import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go

# Function to read the data from CSV files
def read_file(path):
    df = pd.read_csv(path, sep=r'\s+', names=["L", "x", "y"], skiprows=0)
    letters = df["L"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    
    return letters, x, y

# Function to map 'A', 'B', 'C' to colors
def map_labels_to_colors(labels):
    label_colors = {
        'A': 'blue',  # Color for A
        'B': 'green', # Color for B
        'C': 'red'    # Color for C
    }
    return [label_colors[label] for label in labels]

# Folder path where the data files are located
folder_path = 'output/times'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))

# Get the list of frame numbers you want to save
print("Available frames:", [i for i in range(len(files))])
frames_to_save = input("Enter the frame numbers to save (comma-separated, e.g., 0,2,5): ")

# Convert the user input into a list of integers
frames_to_save = [int(i) for i in frames_to_save.split(',')]

# Loop through the specified frames and generate a snapshot (static plot)
for i in frames_to_save:
    if i < len(files):
        letters, x, y = read_file(os.path.join(folder_path, files[i]))
        color_values = map_labels_to_colors(letters)
        
        # Create a static plot for the selected frame (just points, no animation)
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=x, y=y, mode='markers',  # Use 'markers' for just points
            marker=dict(size=8, color=color_values),  # Color based on labels
            name=f"Frame {i}"
        ))
        
        # Update layout for better visualization
        fig.update_layout(
            title=f"Snapshot of Frame {i}",
            xaxis=dict(title="X"),
            yaxis=dict(title="Y"),
            showlegend=False
        )
        
        # Save the figure as an image
        fig.write_image(f"snapshot_frame_{i}.png")
        
        print(f"Snapshot for frame {i} saved as snapshot_frame_{i}.png")
    else:
        print(f"Frame {i} is out of range. Please choose a valid frame number.")

print("Selected snapshots saved!")
