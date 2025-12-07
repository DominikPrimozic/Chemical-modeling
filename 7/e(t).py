# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 17:30:37 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm


directories=["T1","T2","T3","T4","T5","T6","T7","T8"]
tempes={"T1" : 0.25,"T2":0.5,"T3":0.75,"T4":1,"T5":1.25,"T6":1.5,"T7":1.75,"T8":2}
base_path = 'output/MDLJ/'  # Base directory path

# Create a figure for plotting
plt.figure(figsize=(10, 6))

# Define a color map to assign unique colors to each run
cmap = cm.get_cmap("tab10", len(directories))  # Use 'tab10' for a set of distinct colors

# Loop through each directory and plot the data
for idx, directory in enumerate(directories):
    equilibration_file = os.path.join(base_path, directory, 'thermo.txt')
    
    # Load data for equilibration and sampling
    try:
        t,T,K,U,E,p = np.loadtxt(equilibration_file, skiprows=1,unpack=True)
        
        
        # Plot the data for each directory with unique colors
        #plt.plot(t, T, color=cmap(idx), alpha=0.7, linestyle='--')
        plt.plot(t, U, label=f"U {tempes[directory]}", color=cmap(idx), alpha=0.7)
        plt.plot(t, K, linestyle="--",label=f"K {tempes[directory]}", color=cmap(idx), alpha=0.7)
        

    except Exception as e:
        print(f"Error processing files for directory {directory}: {e}")

# Labels and title
plt.xlabel("Time")
plt.ylabel("E")
plt.title("Energy Fluctuation Trends")
plt.legend(title="Temperatures",loc='upper right', bbox_to_anchor=(1, 1), fontsize=9)

plt.tight_layout()  # Make sure everything fits well on the plot
plt.show()
