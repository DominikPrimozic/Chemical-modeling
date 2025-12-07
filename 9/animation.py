# -*- coding: utf-8 -*-
"""
Created on Thu May 22 07:32:23 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# === CONFIGURATION ===
data_folder = '.'        # Folder with the LBM output .dat files
base_filename = 'output/vy'    # Change to 'vx' or 'vy' to visualize those fields
num_frames = 1000        # Number of time steps
Nx, Ny = 100, 50        # Grid size (adjust to match your C++ simulation)

# === LOAD ONE FRAME TO GET SHAPE ===
def load_frame(t):
    filename = os.path.join(data_folder, f"{base_filename}_{t}.dat")
    try:
        data = np.loadtxt(filename)
        return data
    except Exception as e:
        print(f"Failed to load {filename}: {e}")
        return np.zeros((Ny, Nx))  # fallback blank frame

# === SETUP PLOT ===
fig, ax = plt.subplots()
frame0 = load_frame(0)

# Draw the field
im = ax.imshow(frame0, cmap='viridis', origin='lower', interpolation='none')
plt.colorbar(im, ax=ax)
ax.set_title(f'vx at t=0')

# Draw the obstacle as black squares
obstacle_width = 5
half_width = obstacle_width // 2
x_start = Nx//3 - half_width
x_end = Nx//3 + half_width
y_start = 3 * Ny // 8
y_end = 5 * Ny // 8
ax.fill_betweenx(np.arange(y_start, y_end), x_start, x_end, color='black')


# === UPDATE FUNCTION ===
def update(t):
    data = load_frame(t)
    im.set_array(data)
    ax.set_title(f'vx at t={t}')
    return [im]  # Include overlay in returned artists


# === CREATE ANIMATION ===
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=50, blit=False)

# === SHOW OR SAVE ===
# plt.show()  # uncomment to just view it

# === SAVE AS MP4 (requires ffmpeg) ===
#ani.save(f"{base_filename}_animation.mp4", writer='ffmpeg', fps=20)

# === OR SAVE AS GIF (requires imagemagick or pillow) ===
ani.save(f"1000_x3_vy.gif", writer='pillow', fps=20)
