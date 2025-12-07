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
Nx = 100   # Number of grid points in x-direction
Ny = 100   # Number of grid points in y-direction

dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

h = 5.0          # Heat transfer coefficient
phi_fluid = 150  # Ambient temperature

# Dirichlet BC values
left_value = 100.0
right_value = 100.0
bottom_value = 100

# Initialize phi array
phi = np.zeros((Nx, Ny))

# Apply Dirichlet BCs
phi[0, :] = left_value    # Left boundary
phi[-1, :] = right_value  # Right boundary
phi[:, 0] = bottom_value  # Bottom boundary

# Iteration parameters
max_iter = 10000
tolerance = 1e-4
iter = 0
error = 1.0

# Gauss-Seidel iteration with convective BC
while iter < max_iter and error > tolerance:
    error = 0.0
    
    # Update interior points
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            term1 = (phi[i+1, j] + phi[i-1, j]) / dx**2
            term2 = (phi[i, j+1] + phi[i, j-1]) / dy**2
            temp = (term1 + term2) / (2 * (1/dx**2 + 1/dy**2))
            current_error = abs(temp - phi[i, j])
            if current_error > error:
                error = current_error
            phi[i, j] = temp
    
    # Update top boundary (convective BC)
    for i in range(Nx):
        temp = (phi[i, Ny-2] + h * dy * phi_fluid) / (1 + h * dy)
        current_error = abs(temp - phi[i, Ny-1])
        if current_error > error:
            error = current_error
        phi[i, Ny-1] = temp
    
    iter += 1

print(f"Converged in {iter} iterations with error {error:.2e}")

plt.figure(figsize=(10, 8))

img = plt.imshow(phi, cmap="coolwarm", interpolation="bicubic",extent=[0, 1, 0, 1])

# Add colorbar with only min and max ticks
cbar = plt.colorbar(img)
cbar.set_ticks([np.min(phi), np.max(phi)])
cbar.ax.set_yticklabels([f"{np.min(phi):.2f}", f"{np.max(phi):.2f}"])
cbar.ax.set_position([0.8, 0.1, 0.03, 0.8]) 


plt.title("Temperatura preseka", pad=35)

plt.text(0.47, 1.025, f"{phi_fluid}", fontsize=12, ha='left', va='top', transform=plt.gca().transAxes)

# Bottom value annotation (right of the plot)
plt.text(0.47, -0.1, f"{bottom_value}", fontsize=12, ha='left', va='bottom', transform=plt.gca().transAxes)

# Left value annotation (above the plot)
plt.text(-0.05, 0.5, f"{left_value}", fontsize=12, ha='right', va='center', transform=plt.gca().transAxes)

# Right value annotation (above the plot)
plt.text(1.02, 0.5, f"{right_value}", fontsize=12, ha='left', va='center', transform=plt.gca().transAxes)


plt.savefig(f"robin1_{bottom_value:.6f}_{left_value:.6f}.png", dpi=300, bbox_inches='tight')

plt.show()