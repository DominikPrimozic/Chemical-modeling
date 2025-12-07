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
        density = float(f.readline().split('//')[0].strip())
        D = int(f.readline().split('//')[0].strip())
    return N, density, D

# -------------------------
# 🔧 USER INPUT
config_folder = "output/LJ/C/L/"  # Your folder with .bin files
param_file = "input/LJ.txt"
N, rho, dim = read_simulation_params(param_file)
dr = 0.05                    # Bin width

# -------------------------
# 📂 Read all configs from folder
file_paths = sorted([os.path.join(config_folder, f) for f in os.listdir(config_folder) if f.endswith('.bin')])
positions_list = [read_particle_positions(fp, num_particles=N, dimensions=dim) for fp in file_paths]

# -------------------------
# 🚀 Compute g(r)
r, g_r = compute_gr(positions_list, density=rho, dr=dr, dim=dim)

# -------------------------
# 📈 Plot
plt.plot(r, g_r, label="g(r)")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial Distribution Function")
#plt.legend()
plt.tight_layout()
plt.show()
