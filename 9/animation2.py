# -*- coding: utf-8 -*-
"""
Created on Thu May 22 18:57:39 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# === CONFIGURATION ===
data_folder = '.'        # Folder with the LBM output .dat files
num_frames = 1000        # Number of time steps
Nx, Ny = 100, 50         # Grid size (adjust to match your C++ simulation)

# === LOAD VELOCITY MAGNITUDE FRAME ===
def load_magnitude_frame(t):
    try:
        vx = np.loadtxt(os.path.join(data_folder, f"output/vx_{t}.dat"))
        vy = np.loadtxt(os.path.join(data_folder, f"output/vy_{t}.dat"))
        magnitude = np.sqrt(vx**2 + vy**2)
        return magnitude
    except Exception as e:
        print(f"Failed to load data for t={t}: {e}")
        return np.zeros((Ny, Nx))  # fallback blank frame

# === SETUP PLOT ===
fig, ax = plt.subplots()
frame0 = load_magnitude_frame(0)

im = ax.imshow(frame0, cmap='viridis', origin='lower', interpolation='none')
cbar = plt.colorbar(im, ax=ax)
ax.set_title(f'|v| at t=0')

# Draw the obstacle as a filled rectangle
obstacle_width = 5
half_width = obstacle_width // 2
x_start = Nx//3 - half_width
x_end = Nx//3 + half_width
y_start = Ny // 4
y_end = 3 * Ny // 4
y_start = 3 * Ny // 8
y_end = 5 * Ny // 8
# Rectangle patch for obstacle
ax.fill_betweenx(np.arange(y_start, y_end), x_start, x_end, color='black')

# === UPDATE FUNCTION ===
def update(t):
    data = load_magnitude_frame(t)
    im.set_array(data)
    ax.set_title(f'|v| at t={t}')
    return [im]

# === CREATE ANIMATION ===
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=50, blit=False)

# === SHOW OR SAVE ===
# plt.show()

# Save animation
ani.save("Re1000U005_x3_ovira.gif", writer='pillow', fps=20)
