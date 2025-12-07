# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 16:47:44 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 16:14:38 2025

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

def compute_average_velocity_magnitudes(velocities_list):
    """
    Compute the time-averaged magnitudes of the velocities for each particle.

    Args:
        velocities_list: list of (N x dim) arrays of velocities.

    Returns:
        avg_magnitudes: (N,) array of average velocity magnitudes.
    """
    num_snapshots = len(velocities_list)
    N = velocities_list[0].shape[0]  # number of particles

    sum_magnitudes = np.zeros(N)

    for velocities in velocities_list:
        magnitudes = np.linalg.norm(velocities, axis=1)  # (N,)
        sum_magnitudes += magnitudes

    avg_magnitudes = sum_magnitudes / num_snapshots
    return avg_magnitudes

# -------------------------
# 🔧 USER INPUT
velocity_folder = "output/MDLJ/3/vel/"  # Your folder with .bin files (velocity data)
# param_file = "input/LJ.txt"
# N, rho, dim = read_simulation_params(param_file)
N = 100  # Number of particles
rho = 0.4  # Density (not used directly here)
dim = 3  # Dimensions (2D or 3D)
dr = 0.1  # Bin width (not used here for histogram, can be customized)

# -------------------------
# 📂 Read all velocity configurations from folder
file_paths = sorted([os.path.join(velocity_folder, f) for f in os.listdir(velocity_folder) if f.endswith('.bin')])
velocities_list = [read_particle_velocities(fp, num_particles=N, dimensions=dim) for fp in file_paths]

# -------------------------
# 🚀 Compute average velocity magnitudes
avg_velocity_magnitudes = compute_average_velocity_magnitudes(velocities_list)

# -------------------------
# 📈 Plot Histogram of Average Velocity Magnitudes
plt.figure(figsize=(8, 6))
plt.hist(avg_velocity_magnitudes, bins=30, edgecolor='black', color='skyblue', alpha=0.7)
plt.xlabel('Average |v|', fontsize=12)
plt.ylabel('Number of Particles', fontsize=12)
plt.title('Histogram of Time-Averaged Velocity Magnitudes', fontsize=14)
plt.tight_layout()
plt.show()