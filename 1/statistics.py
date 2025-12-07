import os
import numpy as np
import matplotlib.pyplot as plt

# Path to the directory containing your txt files
folder_path = 'output/statistics'

# Initialize lists to store densities and count of 0s
densities = []
count_zeros = []

# Loop through all the files in the directory
for filename in os.listdir(folder_path):
    if filename.endswith('.txt'):
        # Extract the density from the filename (assuming it is the number before .txt)
        
        density = float(filename.replace('.txt', ''))
        densities.append(density)
        
        # Read the file and count the number of 0s
        with open(os.path.join(folder_path, filename), 'r') as file:
            data = file.read().split()  # Read the single line and split by spaces
            count_zeros.append(data.count('0'))  # count '0' entries in the file

# Convert to numpy arrays for easier manipulation
densities = np.array(densities)
count_zeros = np.array(count_zeros)

# Plotting histogram
plt.figure(figsize=(10, 6))
plt.bar(densities, count_zeros, width=0.001, color='blue', edgecolor='black', alpha=0.7)
plt.xlabel('gostota')
plt.ylabel('Uspešni')
plt.title('Število uspešno izvedenih postavitev pri 1000 poskusih')
plt.grid(True)
plt.show()
