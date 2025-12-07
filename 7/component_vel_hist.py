import numpy as np
import matplotlib.pyplot as plt
import os

def read_particle_velocities(filename, num_particles, dimensions=2):
    velocities = np.zeros((num_particles, dimensions))
    with open(filename, "rb") as f:
        for i in range(num_particles):
            velocities[i] = np.fromfile(f, dtype=np.float64, count=dimensions)
    return velocities

# -------------------------
# 🔧 USER INPUT
velocity_folder = "output/MDLJ/3/vel/"  # Your folder with .bin files (velocity data)
N = 100
dim = 3

# -------------------------
# 📂 Read all velocity configurations from folder
file_paths = sorted([os.path.join(velocity_folder, f) for f in os.listdir(velocity_folder) if f.endswith('.bin')])
velocities_list = [read_particle_velocities(fp, num_particles=N, dimensions=dim) for fp in file_paths]

# -------------------------
# 📈 Pick one component (e.g., vx = component 0)
component_index = 0  # 0 for vx, 1 for vy, 2 for vz if 3D

# Extract the chosen component across all configs
all_components = np.concatenate([vel[:, component_index] for vel in velocities_list])

# -------------------------
# 📈 Plot Histogram for the selected velocity component
plt.figure(figsize=(8, 6))
plt.hist(all_components, bins=30, edgecolor='black', color='salmon', alpha=0.7, density=True)
plt.xlabel(f'Velocity Component v[{component_index}]', fontsize=12)
plt.ylabel('Probability Density', fontsize=12)
plt.ylim(0, 1)  # Optional: force y-axis to be between 0 and 1
plt.tight_layout()
plt.show()
