# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 18:52:37 2025

@author: domin
"""


import numpy as np
import matplotlib.pyplot as plt

# Load data
filename = "output/LJ/E/L/equilibration_energy.txt"
_, e1 = np.loadtxt(filename, unpack=True)

filename = "output/LJ/E/L/sampled_energy.txt"
_, e2 = np.loadtxt(filename, unpack=True)

# Combine lengths
total_len = len(e1) + len(e2)

# Normalize x-axis to [0, 1]
x1 = np.linspace(0, len(e1) / total_len, len(e1))
x2 = np.linspace(len(e1) / total_len, 1, len(e2))

# Plot
plt.plot(x1, e1, label="Equilibration", color="blue")
plt.plot(x2, e2, label="Sampling", color="red")

# Labels and title
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Energy Fluctuation")
plt.legend()
plt.show()

# Calculate average energy
avg_E = np.mean(e2)
print(f"Average Energy: {avg_E:.4f}")
print(f"Average Energy per particle: {(avg_E / 1000):.4f}")


