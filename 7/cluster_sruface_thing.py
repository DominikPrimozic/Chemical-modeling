import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def load_positions(file_path, dim=2):
    """ Load particle positions from binary file. """
    with open(file_path, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
        num_particles = data.size // dim
        positions = data.reshape((num_particles, dim))
    return positions

def calculate_pair_distances(positions):
    """ Calculate all unique pair distances. """
    num_particles = positions.shape[0]
    distances = []
    for i in range(num_particles):
        for j in range(i+1, num_particles):
            dist = np.linalg.norm(positions[i] - positions[j])
            distances.append(dist)
    return np.array(distances)

def calculate_potential_energy(positions, epsilon=1.0, sigma=1.0):
    """ Calculate total Lennard-Jones potential energy of the system. """
    num_particles = positions.shape[0]
    potential_energy = 0.0
    for i in range(num_particles):
        for j in range(i+1, num_particles):
            r2 = np.sum((positions[i] - positions[j]) ** 2)
            if r2 < 1e-12:
                r2 = 1e-12
            inv_r6 = (sigma * sigma / r2) ** 3
            potential_energy += 4 * epsilon * (inv_r6 ** 2 - inv_r6)
    return potential_energy

# ==================== Main Script ====================

# --- User Settings ---
pos_folder = "output/cluster/5/pos/"     # Folder containing position .bin files
dim = 2                 # Dimension of the system (2 or 3)
num_bins = 50           # Number of bins for distance histogram
max_distance = 10.0     # Maximum distance for histogram (adjust if needed)

# --- Load all position files ---
position_files = sorted(glob.glob(os.path.join(pos_folder, "*.bin")))
print(f"Found {len(position_files)} position files.")

energies = []
distance_means = []

# --- Process each position file ---
for file_path in position_files:
    positions = load_positions(file_path, dim=dim)
    
    distances = calculate_pair_distances(positions)
    potential_energy = calculate_potential_energy(positions)
    
    energies.append(potential_energy)
    #distance_means.append(np.mean(distances))  # Mean pair distance
    distance_means.append(np.median(distances))

# --- Convert to arrays ---
energies = np.array(energies)
distance_means = np.array(distance_means)

# --- Plot Energy vs Mean Pair Distance ---
plt.figure(figsize=(8,6))
plt.scatter(distance_means, energies, c='red', edgecolors='k',s=10)
plt.xlabel('Mean Pair Distance')
plt.ylabel('Potential Energy')
plt.title('Potential Energy vs Mean Pair Distance (Free Cluster)')
plt.tight_layout()
plt.show()


import plotly.graph_objects as go
import plotly.io as pio

# --- Your function (already correct) ---
def create_circle(cx, cy, r=0.5, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return x, y

# --- Plotly setup ---
fig = go.Figure()

def add_circles_for_configuration(coords):
    traces = []
    for cx, cy in coords:
        x, y = create_circle(cx, cy, r=d/2)
        traces.append(go.Scatter(
            x=x, y=y,
            fill="none",  
            mode="lines",
            line=dict(color='blue', dash='dot'),
            opacity=0.6,
            showlegend=False
        ))
        traces.append(go.Scatter(
            x=[cx], y=[cy],
            mode="markers",
            marker=dict(size=2, color="red", symbol="circle"),
            showlegend=False
        ))
    return traces

# --- Find lowest energy configuration ---
min_index = np.argmin(energies)
lowest_energy_file = position_files[min_index]
lowest_positions = load_positions(lowest_energy_file, dim=dim)

# --- Now use your plotter on lowest_positions ---
# (lowest_positions should have shape (N, 2) if 2D)

# Assume particle diameter d
d = 1.0  # or whatever your particle size is

# Add all circles
for trace in add_circles_for_configuration(lowest_positions):
    fig.add_trace(trace)

# Set plot layout nicely
fig.update_layout(
    title=f"Lowest Energy Configuration (E={energies[min_index]:.3f})",
    xaxis_title="x",
    yaxis_title="y",
    yaxis_scaleanchor="x",  # equal aspect ratio
    template="plotly_white",
    width=600,
    height=600
)

pio.show(fig, renderer='browser')
