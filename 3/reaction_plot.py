# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 20:27:45 2025

@author: domin
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt

wanted=0

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

profile=f"output/reaction/profile/{values['D']:.6f}_{values['a']:.6f}_{values['c0']:.6f}_{values['T0']:.6f}_{values['TH']:.6f}.txt"
flux=f"output/reaction/flux/{values['D']:.6f}_{values['a']:.6f}_{values['c0']:.6f}_{values['T0']:.6f}_{values['q']:.6f}.txt"
isolated=f"output/reaction/isolated/{values['D']:.6f}_{values['k']:.6f}_{values['a']:.6f}_{values['H']:.6f}_{values['c0']:.6f}.txt"

if wanted==0:
    subprocess.run(["run_reactor_isolated.exe"], check=True)
    x,c,dc=np.loadtxt(isolated,unpack=True)
    
    plt.figure(figsize=(8, 5))
    plt.plot(x, c, linestyle="-", color="royalblue", label="c")
    plt.xlabel("x")
    plt.ylabel("c")
    plt.title("Isolated")
    plt.legend()
    plt.show()

elif wanted==1:
    subprocess.run(["run_reactor_flux.exe"], check=True)
    x,c,dc,T,dT=np.loadtxt(flux,unpack=True)
    
    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True)
    
    axs[0].plot(x, c, linestyle="-", color="darkorange", label="c")
    axs[0].set_ylabel("c")
    axs[0].set_title("Flux Concentration")
    axs[0].legend()
    
    
    axs[1].plot(x, T, linestyle="-", color="crimson", label="T")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("T")
    axs[1].set_title("Flux Temperature")
    axs[1].legend()
    
    
    plt.tight_layout()
    plt.show()

else: 
    subprocess.run(["run_reactor_profile.exe"], check=True)
    x,c,dc=np.loadtxt(profile,unpack=True)
    
    plt.figure(figsize=(8, 5))
    plt.plot(x, c, linestyle="-", color="seagreen", label="c")
    plt.xlabel("x")
    plt.ylabel("c")
    plt.title("Profile")
    plt.legend()
    
    plt.show()

#