# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 13:02:19 2025

@author: domin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.offline as pyo
import subprocess
import seaborn as sns
def read_input_file(filename):
    values = {}
    keys = [
        "left", 
        "right", 
        "bottom", 
        "top", 
        "source"
    ]

    with open(filename, 'r') as f:
        lines = f.readlines()
        
        for key, line in zip(keys, lines):
            number = int(line.split("//")[0].strip())
            values[key] = number

    return values

# Example usage
input_data = read_input_file("input/poisson/input.txt")
result = subprocess.run(["block_poisson.exe"], capture_output=True, text=True)

path=f"output/poisson/{input_data['bottom']:.6f}_{input_data['left']:.6f}.txt"

data=np.loadtxt(path, delimiter="\t")

rows, cols = data.shape
new_data = np.zeros((rows + 2, cols + 2))  # Create a larger matrix

# Fill the inner part with original data
new_data[1:-1, 1:-1] = data  

# Assign corner values
new_data[-1, :] = input_data['top']  # Top row with x0
new_data[0, :] = input_data['bottom']  # Bottom row with xn
new_data[:, 0] = input_data['left']  # Left column with y0
new_data[:, -1] = input_data['right']
new_data=np.flipud(new_data)
plt.figure(figsize=(10, 8))

img = plt.imshow(new_data, cmap="coolwarm", interpolation="bicubic",extent=[0, 1, 1, 0])

# Add colorbar with only min and max ticks
cbar = plt.colorbar(img)
cbar.set_ticks([np.min(new_data), np.max(new_data)])
cbar.ax.set_yticklabels([f"{np.min(new_data):.2f}", f"{np.max(new_data):.2f}"])
cbar.ax.set_position([0.8, 0.1, 0.03, 0.8]) 


plt.title(f"Temperatura preseka (S={input_data['source']})", pad=35)

plt.text(0.47, 1.025, f"{input_data['top']}", fontsize=12, ha='left', va='top', transform=plt.gca().transAxes)

# Bottom value annotation (right of the plot)
plt.text(0.47, -0.1, f"{input_data['bottom']}", fontsize=12, ha='left', va='bottom', transform=plt.gca().transAxes)

# Left value annotation (above the plot)
plt.text(-0.05, 0.5, f"{input_data['left']}", fontsize=12, ha='right', va='center', transform=plt.gca().transAxes)

# Right value annotation (above the plot)
plt.text(1.02, 0.5, f"{input_data['right']}", fontsize=12, ha='left', va='center', transform=plt.gca().transAxes)


plt.savefig(f"laplace_{input_data['bottom']:.6f}_{input_data['source']:.6f}.png", dpi=300, bbox_inches='tight')

plt.show()