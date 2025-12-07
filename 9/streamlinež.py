# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:08:15 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Streamline plot for velocity field at last frame
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# === CONFIGURATION ===
data_folder = '.'        # Folder with the LBM output .dat files
Nx, Ny = 100, 50         # Grid size (match simulation)
frame = 999              # Last frame if num_frames = 1000

# === LOAD VX AND VY ===
try:
    vx = np.loadtxt(os.path.join(data_folder, f"output/vx_{frame}.dat"))
    vy = np.loadtxt(os.path.join(data_folder, f"output/vy_{frame}.dat"))
except Exception as e:
    print(f"Failed to load velocity data: {e}")
    vx = np.zeros((Ny, Nx))
    vy = np.zeros((Ny, Nx))

# === SETUP GRID ===
x = np.linspace(0, Nx-1, Nx)
y = np.linspace(0, Ny-1, Ny)
X, Y = np.meshgrid(x, y)

# === CREATE PLOT ===
fig, ax = plt.subplots(figsize=(10, 5))
strm = ax.streamplot(X, Y, vx, vy, color=np.sqrt(vx**2 + vy**2), cmap='plasma', linewidth=1)

# Add colorbar
cbar = plt.colorbar(strm.lines)
cbar.set_label('|v|')

# === DRAW OBSTACLE RECTANGLE ===
obstacle_width = 5
half_width = obstacle_width // 2
x_start = Nx//3 - half_width
x_end = Nx//3 + half_width
y_start = 3 * Ny // 8
y_end = 5 * Ny // 8

# Draw filled rectangle for obstacle
ax.fill_betweenx(np.arange(y_start, y_end), x_start, x_end, color='black')

# === FINALIZE PLOT ===
ax.set_title(f'Streamlines at t={frame}')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim(0, Nx-1)
ax.set_ylim(0, Ny-1)
ax.set_aspect('equal')

plt.tight_layout()
plt.show()