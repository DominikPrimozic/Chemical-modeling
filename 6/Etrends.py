import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm

# Directories for each dataset
directories = ['1', '2', '3', '4', '5']
dense={'1':0.1, '2':0.2, '3':0.3, '4':0.4, '5':0.5,}
directories=["05","10","20","50","100"]
tempes={"05":0.5,"10":1,"2":2,"50":5,"100":10}
base_path = 'output/LJ/E/T/'  # Base directory path

# Create a figure for plotting
plt.figure(figsize=(10, 6))

# Define a color map to assign unique colors to each run
cmap = cm.get_cmap("tab10", len(directories))  # Use 'tab10' for a set of distinct colors

# Loop through each directory and plot the data
for idx, directory in enumerate(directories):
    equilibration_file = os.path.join(base_path, directory, 'equilibration_energy.txt')
    sampled_file = os.path.join(base_path, directory, 'sampled_energy.txt')
    
    # Load data for equilibration and sampling
    try:
        _, e1 = np.loadtxt(equilibration_file, unpack=True)
        _, e2 = np.loadtxt(sampled_file, unpack=True)
        
        # Combine lengths
        total_len = len(e1) + len(e2)
        
        # Normalize x-axis to [0, 1]
        x1 = np.linspace(0, len(e1) / total_len, len(e1))
        x2 = np.linspace(len(e1) / total_len, 1, len(e2))
        
        # Plot the data for each directory with unique colors
        plt.plot(x1, e1, color=cmap(idx), alpha=0.7, linestyle='--')
        #plt.plot(x2, e2, label=f"{dense[directory]}", color=cmap(idx), alpha=0.7)
        plt.plot(x2, e2, label=f"{tempes[directory]}", color=cmap(idx), alpha=0.7)
        
        # Add a vertical line at the end of equilibration (x1[-1] is where equilibration ends)
        plt.axvline(x=x1[-1], color=cmap(idx), linestyle=':', linewidth=2)

    except Exception as e:
        print(f"Error processing files for directory {directory}: {e}")

# Labels and title
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Energy Fluctuation Trends")
plt.legend(title="Temperatures",loc='upper right', bbox_to_anchor=(1, 1), fontsize=9)

plt.tight_layout()  # Make sure everything fits well on the plot
plt.show()
