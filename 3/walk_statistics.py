# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 11:12:47 2025

@author: domin
"""
import numpy as np
import matplotlib.pyplot as plt

N=100
d=0.300000

"""
file = f"output\\random_walk\\steps\\{N}_{d:.6f}_{1}.txt"
with open(file, "r") as f:
    average_failed = float(f.readline().strip())

for rep in range(1,1000):
    file = f"output\\random_walk\\steps\\{N}_{d:.6f}_{rep}.txt"
    failed_steps=np.loadtxt(file,unpack=True, skiprows=2)
    average_failed+=failed_steps
    
average_failed/=998
    


file = f"output\\random_walk\\distance\\{N}_{d:.6f}_{1}.txt"
average_distance = np.loadtxt(file,unpack=True)
    
for rep in range(2,1000):
    file = f"output\\random_walk\\distance\\{N}_{d:.6f}_{rep}.txt"
    failed_steps=np.loadtxt(file,unpack=True)
    average_distance+=failed_steps
    
average_distance/=998
    

plt.plot(average_distance)
plt.show()

distances_array=[]

for i in range(1,100):
    d=i*0.006
    file = f"output\\random_walk\\distance\\{N}_{d:.6f}_{1}.txt"
    average_distanced = np.loadtxt(file,unpack=True)
    for rep in range(2,500):
        file = f"output\\random_walk\\distance\\{N}_{d:.6f}_{rep}.txt"
        failed_steps=np.loadtxt(file,unpack=True)
        average_distanced+=failed_steps
    average_distanced/=498
    distances_array.append(average_distanced)
    
distances_array=np.array(distances_array)
"""

"""
first part
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import seaborn as sns

# Initialize storage
N = 100  # Adjust as needed
distances_array = []
d_values = []

# Compute distances
for i in range(1, 100):
    d = i * 0.006
    d_values.append(d)
    file = f"output/random_walk/distance/{N}_{d:.6f}_1.txt"
    average_distanced = np.loadtxt(file, unpack=True)
    
    for rep in range(2, 300):
        file = f"output/random_walk/distance/{N}_{d:.6f}_{rep}.txt"
        failed_steps = np.loadtxt(file, unpack=True)
        average_distanced += failed_steps
    
    average_distanced /= 398  # Averaging over repetitions
    distances_array.append(average_distanced)

distances_array = np.array(distances_array)

# Reverse the density axis
d_values = d_values[::-1]
distances_array = distances_array[::-1]

# Create Potential Energy Surface Plot (3D)
fig_pes = go.Figure()
x = np.arange(distances_array.shape[1])  # Assuming steps on x-axis
y = np.array(d_values)  # Density steps
X, Y = np.meshgrid(x, y)

fig_pes.add_trace(go.Surface(z=distances_array, x=X, y=Y, colorscale='Viridis'))
fig_pes.update_layout(title='Oddaljenost', scene=dict(xaxis_title='Step', yaxis_title='Density', zaxis_title='Oddaljenost od izhodišča'))
pio.show(fig_pes, renderer="browser")

# Create Heatmap
plt.figure(figsize=(10, 6))
ax = sns.heatmap(distances_array, cmap='viridis', xticklabels=50, yticklabels=10)
plt.xlabel("Step Index")
plt.ylabel("Density")
plt.title("Potential Energy Heatmap")
plt.xlim(0, max(x))
plt.ylim(0, max(y))

# Fix density tick labels
ytick_positions = np.linspace(0, len(y) - 1, num=5, dtype=int)  # Pick 10 evenly spaced indices
ytick_labels = [f"{y[i]:.3f}" for i in ytick_positions[::-1]]  # Format density values
plt.yticks(ytick_positions, ytick_labels)

plt.show()
"""

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import seaborn as sns

# Initialize storage
N = 100  # Adjust as needed
failed_array = []
d_values = []
pio.renderers.default = "browser"
# Compute distances
for i in range(1, 50):
    d = i * 0.012
    d_values.append(d)
    failed_step=[]
    for j in range(1,50):
        step=j*0.04
        file = f"output/random_walk/steps/{N}_{d:.6f}_{step:.6f}_0.txt"
        average_failed_step = np.loadtxt(file, unpack=True, skiprows=2) #250 total steps
    
        for rep in range(1, 100):
            file = f"output/random_walk/steps/{N}_{d:.6f}_{step:.6f}_{rep}.txt"
            failed_steps = np.loadtxt(file, unpack=True, skiprows=2)
            average_failed_step += failed_steps
        average_failed_step /= 100  
        failed_step.append(average_failed_step)
    failed_array.append(failed_step)

failed_array = np.array(failed_array)
d_values = np.array(d_values)
step_values = np.arange(1, 50) * 0.04  # Step sizes

# Select the index corresponding to step=50 (since index starts from 0, it's at 49)
step_idx = 49
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

# --- Heatmap ---
plt.figure(figsize=(10, 8))
ax = sns.heatmap(failed_slice, cmap="viridis", cbar_kws={'label': 'Failed Attempts'})

# Reduce number of ticks
num_ticks = 5
x_indices = np.linspace(0, len(d_values) - 1, num_ticks, dtype=int)
y_indices = np.linspace(0, len(step_values) - 1, num_ticks, dtype=int)

ax.set_xticks(x_indices)
ax.set_xticklabels(np.round(d_values[x_indices], 3))

ax.set_yticks(y_indices)
ax.set_yticklabels(np.round(step_values[y_indices], 3))

# Labels & Title
plt.xlabel("Density")
plt.ylabel("Step Size")
plt.title("Heatmap of Failed Attempts")

plt.show()

#######################################################
step_idx = 1
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 1)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 10
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 10)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 20
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 20)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 30
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 30)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 40
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 40)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 50
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 50)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 100
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 100)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser

step_idx = 200
failed_slice = failed_array[:, :, step_idx]  # Shape: (49, 49)

# --- 3D Surface Plot ---
fig = go.Figure(data=[go.Surface(z=failed_slice, x=d_values, y=step_values)])
fig.update_layout(title="3D Surface Plot (Density vs Step Size vs Failed Attempts at 200)",
                  scene=dict(xaxis_title="Density", yaxis_title="Step Size", zaxis_title="Failed Attempts"),
                  autosize=True)
fig.show()  # This will open the plot in your browser