# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 16:14:38 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 20:11:40 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def read_particle_velocities(filename, num_particles, dimensions=2):
    """
    Read a single configuration of particle velocities from a binary file.
    Returns: (N x dim) array of velocities.
    """
    velocities = np.zeros((num_particles, dimensions))
    with open(filename, "rb") as f:
        for i in range(num_particles):
            velocities[i] = np.fromfile(f, dtype=np.float64, count=dimensions)
    return velocities

def compute_velocity_magnitudes(velocities_list):
    """
    Compute the magnitudes of the velocities for each configuration.

    Args:
        velocities_list: list of (N x dim) arrays of velocities.

    Returns:
        magnitudes: List of velocity magnitudes (N, )
    """
    magnitudes_list = []

    for velocities in velocities_list:
        magnitudes = np.linalg.norm(velocities, axis=1)  # Compute the magnitude of each velocity
        magnitudes_list.append(magnitudes)
    
    return magnitudes_list

# -------------------------
# 🔧 USER INPUT
velocity_folder = "output/MDLJ/3/vel/"  # Your folder with .bin files (velocity data)
#param_file = "input/LJ.txt"
#N, rho, dim = read_simulation_params(param_file)
N = 100  # Number of particles
rho = 0.4  # Density (not used directly here)
dim = 3  # Dimensions (2D or 3D)
dr = 0.1  # Bin width (not used here for histogram, can be customized)

# -------------------------
# 📂 Read all velocity configurations from folder
file_paths = sorted([os.path.join(velocity_folder, f) for f in os.listdir(velocity_folder) if f.endswith('.bin')])
velocities_list = [read_particle_velocities(fp, num_particles=N, dimensions=dim) for fp in file_paths]

# -------------------------
# 🚀 Compute velocity magnitudes
velocity_magnitudes_over_time = compute_velocity_magnitudes(velocities_list)

# Flatten the list of velocity magnitudes for histogram
all_velocity_magnitudes = np.concatenate(velocity_magnitudes_over_time)

# -------------------------
# 📈 Plot Histogram of Velocity Magnitudes
# 📈 Plot Histogram of Velocity Magnitudes
plt.figure(figsize=(8, 6))
plt.hist(all_velocity_magnitudes, bins=30, edgecolor='black', color='skyblue', alpha=0.7, density=True)
plt.xlabel('|v|', fontsize=12)
plt.ylabel('Probability Density', fontsize=12)
plt.ylim(0, 1)  # Optional: Force y-axis from 0 to 1 if you want
plt.tight_layout()
plt.show()
