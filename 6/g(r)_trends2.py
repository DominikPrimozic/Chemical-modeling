import numpy as np
import matplotlib.pyplot as plt
import os

def read_particle_positions(filename, num_particles, dimensions=2):
    """
    Read a single configuration from a binary file.
    Returns: (N x dim) array
    """
    positions = np.zeros((num_particles, dimensions))
    with open(filename, "rb") as f:
        for i in range(num_particles):
            positions[i] = np.fromfile(f, dtype=np.float64, count=dimensions)
    return positions

def compute_gr(positions_list, density, dr=0.01, r_max=None, dim=2):
    """
    Compute radial distribution function g(r).

    Args:
        positions_list: list of (N x dim) arrays
        density: number density (N / volume)
        dr: bin width
        r_max: max radius to compute (default: half box)
        dim: 1, 2, or 3

    Returns:
        r: array of r values
        g_r: g(r) values
    """
    num_configs = len(positions_list)
    N = positions_list[0].shape[0]
    volume = N / density
    box_length = volume ** (1 / dim)

    if r_max is None:
        r_max = box_length / 2.0

    N_bins = int(r_max / dr)
    hist = np.zeros(N_bins)

    for pos in positions_list:
        for i in range(N):
            for j in range(i + 1, N):
                rij = pos[i] - pos[j]
                rij -= box_length * np.round(rij / box_length)  # minimum image
                r = np.linalg.norm(rij)
                bin_index = int(r / dr)
                if bin_index < N_bins:
                    hist[bin_index] += 2

    r_vals = (np.arange(N_bins) + 0.5) * dr
    if dim == 3:
        shell_volume = 4 * np.pi * r_vals**2 * dr
    elif dim == 2:
        shell_volume = 2 * np.pi * r_vals * dr
    elif dim == 1:
        shell_volume = 2 * dr
    else:
        raise ValueError("Unsupported dimension")

    g_r = hist / (num_configs * N * density * shell_volume)
    return r_vals, g_r

def read_simulation_params(filepath):
    with open(filepath, 'r') as f:
        N = int(f.readline().split('//')[0].strip())
        D = int(f.readline().split('//')[0].strip())
    return N, D

# -------------------------
# 🔧 USER INPUT
base_config_folder = "output/LJ/C/T/"  # Base folder with subdirectories for different runs
param_file = "input/LJ.txt"
N=100
dim=2
dr = 0.1                     # Bin width
rho=0.4


# -------------------------
# 📂 Read all directories and compute g(r) for each
directories=["05","10","20","50","100"]
tempes={"05":0.5,"10":1,"20":2,"50":5,"100":10}
plt.figure(figsize=(8, 6))  # Plot settings

for directory in directories:
    config_folder = os.path.join(base_config_folder, directory)
    file_paths = sorted([os.path.join(config_folder, f) for f in os.listdir(config_folder) if f.endswith('.bin')])
    positions_list = [read_particle_positions(fp, num_particles=N, dimensions=dim) for fp in file_paths]
    
    # Get the corresponding density for this directory
    
    # Compute g(r) for the current directory
    r, g_r = compute_gr(positions_list, density=rho, dr=dr, dim=dim)
    
    # Plot g(r) for this directory
    plt.plot(r, g_r, label=f"(T={tempes[directory]})")

# -------------------------
# 📈 Final Plot
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial Distribution Function")
plt.legend(title="Runs", loc="best")
plt.xlim(0,5)
plt.tight_layout()
plt.show()
