# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 18:55:19 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
Lx = 1.0
Ly = 1.0
Nx = 100
Ny = 100
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

h = 5.0             # Heat transfer coefficient
phi_fluid = 150     # Ambient temperature
S = 1000.0           # Uniform source term (W/m^3 or similar unit)

# Dirichlet BCs
left_value = 50
bottom_value = 100
top_value = 100.0  # top = j = Ny - 1


# Initialize phi
phi = np.zeros((Nx, Ny))

# Apply Dirichlet BCs
phi[0, :] = left_value
phi[:, 0] = bottom_value
phi[:, -1] = top_value

# Iteration params
max_iter = 10000
tolerance = 1e-4
iter = 0
error = 1.0

# Coefficients
inv_dx2 = 1.0 / dx**2
inv_dy2 = 1.0 / dy**2
denom = 2 * (inv_dx2 + inv_dy2)

# Gauss-Seidel loop
while iter < max_iter and error > tolerance:
    error = 0.0

    # Update interior points with source term
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            term1 = (phi[i + 1, j] + phi[i - 1, j]) * inv_dx2
            term2 = (phi[i, j + 1] + phi[i, j - 1]) * inv_dy2
            temp = (term1 + term2 + S) / denom
            current_error = abs(temp - phi[i, j])
            if current_error > error:
                error = current_error
            phi[i, j] = temp

    # Apply convective BC on right (i = Nx - 1)
    for j in range(Ny):
        temp = (phi[Nx - 2, j] + h * dx * phi_fluid) / (1 + h * dx)
        current_error = abs(temp - phi[Nx - 1, j])
        if current_error > error:
            error = current_error
        phi[Nx - 1, j] = temp

    iter += 1

print(f"Converged in {iter} iterations with error {error:.2e}")

# Plot
plt.figure(figsize=(10, 8))
img = plt.imshow(phi, cmap="coolwarm", interpolation="bicubic", extent=[0, 1, 0, 1], origin='lower')

# Colorbar
cbar = plt.colorbar(img)
cbar.set_ticks([np.min(phi), np.max(phi)])
cbar.ax.set_yticklabels([f"{np.min(phi):.2f}", f"{np.max(phi):.2f}"])
cbar.ax.set_position([0.8, 0.1, 0.03, 0.8])

plt.title("Temperatura preseka", pad=35)

# Boundary annotations
plt.text(0.47, -0.1, f"{left_value}", fontsize=12, ha='left', va='bottom', transform=plt.gca().transAxes)
plt.text(0.47, 1.025, f"{phi_fluid}", fontsize=12, ha='left', va='top', transform=plt.gca().transAxes)
plt.text(-0.05, 0.5, f"{bottom_value}", fontsize=12, ha='right', va='center', transform=plt.gca().transAxes)
plt.text(1.02, 0.5, f"{top_value}", fontsize=12, ha='left', va='center', transform=plt.gca().transAxes)

plt.savefig(f"robin_source_{S:.1f}.png", dpi=300, bbox_inches='tight')
plt.show()
