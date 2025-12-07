import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go

def read_file(path):
    df = pd.read_csv(path, sep=r'\s+', names=["L", "x", "y"], skiprows=0)
    return df["L"].to_numpy(), df["x"].to_numpy(), df["y"].to_numpy()

def map_labels_to_colors(labels):
    label_colors = {'A': 'blue', 'B': 'green', 'C': 'red', 'P': 'gray'}
    return [label_colors[label] for label in labels]

# Folder path and sorting files
folder_path = 'output/porous/times'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))

# Read first frame
letters, x, y = read_file(os.path.join(folder_path, files[0]))

# Extract P (static)
p_indices = np.where(letters == 'P')[0]
p_x, p_y = x[p_indices], y[p_indices]

# Extract A, B, C (moving)
moving_indices = np.where(letters != 'P')[0]
moving_x, moving_y = x[moving_indices], y[moving_indices]
moving_colors = map_labels_to_colors(letters[moving_indices])

# Create figure
fig = go.Figure()

# **Static P points**
if len(p_x) > 0:
    fig.add_trace(go.Scatter(
        x=p_x, y=p_y,
        mode='markers',
        marker=dict(size=8, color='gray'),
        name='P (Static)'
    ))

# **Add first frame of A, B, C (which will animate)**
scatter_trace = go.Scatter(
    x=moving_x, y=moving_y,
    mode='markers',
    marker=dict(size=8, color=moving_colors),
    name='Particles (A, B, C)'
)
fig.add_trace(scatter_trace)

# Create animation frames
frames = []
for file_name in files:
    letters, x, y = read_file(os.path.join(folder_path, file_name))

    moving_indices = np.where(letters != 'P')[0]  # Exclude P
    frame_x, frame_y = x[moving_indices], y[moving_indices]
    frame_colors = map_labels_to_colors(letters[moving_indices])

    frames.append(go.Frame(data=[go.Scatter(
        x=frame_x, y=frame_y,
        mode='markers',
        marker=dict(size=8, color=frame_colors),
    )]))

fig.frames = frames

# Add Play/Pause Buttons and Slider
fig.update_layout(
    title="Animated Reaction (P Static, A/B/C Moving)",
    xaxis=dict(title="X", range=[-10, 10]),
    yaxis=dict(title="Y", range=[-10, 10]),
    showlegend=True,
    updatemenus=[dict(
        type="buttons",
        showactive=True,
        buttons=[
            dict(label="Play", method="animate", args=[None, dict(frame=dict(duration=200, redraw=True), fromcurrent=True)]),
            dict(label="Pause", method="animate", args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])
        ]
    )],
    sliders=[{
        "active": 0,
        "steps": [{"args": [[str(i)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                   "label": str(i),
                   "method": "animate"} for i in range(len(files))],
        "currentvalue": {"prefix": "Frame: ", "visible": True}
    }]
)

fig.show()
fig.write_html("reaction_animation_static_P.html")
