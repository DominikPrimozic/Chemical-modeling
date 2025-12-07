# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 20:27:45 2025

@author: domin
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt

wanted=0

def modify_line(filename, line_number, new_value):
    with open(filename, 'r') as file:
        lines = file.readlines()

    if 1 <= line_number <= len(lines):
        parts = lines[line_number - 1].split("//")
        lines[line_number - 1] = f"{new_value} //{parts[1]}" if len(parts) > 1 else f"{new_value}\n"

        with open(filename, 'w') as file:
            file.writelines(lines)
        print(f"Line {line_number} updated successfully.")
    else:
        print("Error: Line number out of range.")

def read_input_file(filename):
    data = {}
    variable_names = [
        "k", "c0", "D", "H", "a", "l", "T0", "TH", "q",
        "Ea", "Ed", "heat", "dx", "h", "s0"
    ]

    with open(filename, "r") as file:
        index = 0
        for line in file:
            line = line.split("//")[0].strip()  # Remove comments and whitespace
            if line:  # Only process non-empty lines
                try:
                    data[variable_names[index]] = float(line)
                    index += 1
                except ValueError:
                    print(f"Warning: Could not convert '{line}' to a number.")
                if index >= len(variable_names):  # Stop if we have all variables
                    break

    return data

input_file = "input/input.txt"
values = read_input_file(input_file)


def plot_trends(wanted, values, label=None, color=None):
    profile = f"output/reaction/profile/{values['D']:.6f}_{values['a']:.6f}_{values['c0']:.6f}_{values['T0']:.6f}_{values['TH']:.6f}.txt"
    flux = f"output/reaction/flux/{values['D']:.6f}_{values['a']:.6f}_{values['c0']:.6f}_{values['T0']:.6f}_{values['q']:.6f}.txt"
    isolated = f"output/reaction/isolated/{values['D']:.6f}_{values['k']:.6f}_{values['a']:.6f}_{values['H']:.6f}_{values['c0']:.6f}.txt"


    default_color = "royalblue"

    # Check if color is valid (i.e., not None, not an array, etc.)
    if isinstance(color, (list, np.ndarray)):
        plot_color = color
    else:
        plot_color = color if color is not None else default_color
        
        
    if wanted == 0:
        subprocess.run(["run_reactor_isolated.exe"], check=True)
        x, c, dc = np.loadtxt(isolated, unpack=True)
        
        ax = plt.gca()  # Get current axis
        ax.plot(x, c, linestyle="-", color=plot_color, label=label or "c")
        ax.set_xlabel("x")
        ax.set_ylabel("c")
        ax.set_title("Isolated")
        ax.legend()

    elif wanted == 1:
        subprocess.run(["run_reactor_flux.exe"], check=True)
        x, c, dc, T, dT = np.loadtxt(flux, unpack=True)
        
        fig = plt.gcf()  # Get current figure (creates one if needed)
        if not fig.get_axes():  # If no axes exist, create them
            fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True)
        else:
            axs = fig.get_axes()  # Retrieve existing axes

        axs[0].plot(x, c, linestyle="-", color=plot_color, label=label or "c")
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("c")
        axs[0].set_title("Flux Concentration")
        axs[0].legend()
        
        axs[1].plot(x, T, linestyle="-", color=plot_color, label=label or "T")
        axs[1].set_xlabel("x")
        axs[1].set_ylabel("T")
        axs[1].set_title("Flux Temperature")
        axs[1].legend()

    else: 
        subprocess.run(["run_reactor_profile.exe"], check=True)
        x, c, dc = np.loadtxt(profile, unpack=True)
        
        ax = plt.gca()  # Get current axis
        ax.plot(x, c, linestyle="-", color=plot_color, label=label or "c")
        ax.set_xlabel("x")
        ax.set_ylabel("c")
        ax.set_title("Profile")
        ax.legend()

plt.figure(figsize=(8, 5))        
#fig, axs = plt.subplots(1, 2, figsize=(10, 5))

values = read_input_file(input_file)
#plot_trends(1, values)

colors = plt.cm.viridis(np.linspace(0, 1, 11))
for i in range(0,11):
    r=i*0.3
    modify_line(input_file, 5, r)
    values = read_input_file(input_file)
    col = colors[i] 
    plot_trends(2, values, label=f"{r:.1f}",color=col)
    
plt.savefig(f"profile_a_H10.png")