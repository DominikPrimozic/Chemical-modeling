# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 16:21:05 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
Lx = 1.0  # Length in x-direction
Ly = 1.0  # Length in y-direction
Nx = 100  # Number of grid points in x-direction
Ny = 100  # Number of grid points in y-direction

dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

h = 1.0          # Heat transfer coefficient
phi_fluid = 150  # Ambient temperature

# Dirichlet BC values
left_value = 50
bottom_value = 100
top_value = 100.0  # top = j = Ny - 1

# Initialize phi array
phi = np.zeros((Nx, Ny))

# Apply Dirichlet BCs
phi[0, :] = left_value       # Left boundary (i=0)
phi[:, 0] = bottom_value     # Bottom boundary (j=0)
phi[:, -1] = top_value       # Top boundary (j=Ny-1)

# Iteration parameters
max_iter = 10000
tolerance = 1e-4
iter = 0
error = 1.0

# Gauss-Seidel iteration with convective BC on the RIGHT
while iter < max_iter and error > tolerance:
    error = 0.0

    # Update interior points
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            term1 = (phi[i + 1, j] + phi[i - 1, j]) / dx**2
            term2 = (phi[i, j + 1] + phi[i, j - 1]) / dy**2
            temp = (term1 + term2) / (2 * (1 / dx**2 + 1 / dy**2))
            current_error = abs(temp - phi[i, j])
            if current_error > error:
                error = current_error
            phi[i, j] = temp

    # Apply convective BC on RIGHT edge (i = Nx - 1)
    for j in range(Ny):
        temp = (phi[Nx - 2, j] + h * dx * phi_fluid) / (1 + h * dx)
        current_error = abs(temp - phi[Nx - 1, j])
        if current_error > error:
            error = current_error
        phi[Nx - 1, j] = temp

    iter += 1

print(f"Converged in {iter} iterations with error {error:.2e}")

# Plotting
plt.figure(figsize=(10, 8))

img = plt.imshow(phi, cmap="coolwarm", interpolation="bicubic", extent=[0, 1, 0, 1], origin='lower')

# Colorbar
cbar = plt.colorbar(img)
cbar.set_ticks([np.min(phi), np.max(phi)])
cbar.ax.set_yticklabels([f"{np.min(phi):.2f}", f"{np.max(phi):.2f}"])
cbar.ax.set_position([0.8, 0.1, 0.03, 0.8])

plt.title("Temperatura preseka", pad=35)

# Annotations (updated to reflect boundaries)
plt.text(0.47, -0.1, f"{left_value}", fontsize=12, ha='left', va='bottom', transform=plt.gca().transAxes)
plt.text(0.47, 1.025, f"{phi_fluid}", fontsize=12, ha='left', va='top', transform=plt.gca().transAxes)
plt.text(-0.05, 0.5, f"{bottom_value}", fontsize=12, ha='right', va='center', transform=plt.gca().transAxes)
plt.text(1.02, 0.5, f"{top_value}", fontsize=12, ha='left', va='center', transform=plt.gca().transAxes)

plt.savefig(f"robin_1_{bottom_value:.6f}_{left_value:.6f}.png", dpi=300, bbox_inches='tight')
plt.show()
