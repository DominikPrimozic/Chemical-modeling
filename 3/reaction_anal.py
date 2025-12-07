# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 09:03:35 2025

@author: domin
"""

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

def zero(x,values):
    return values["c0"]+(values["k"]/values["D"])/2 * (x*x - values["H"]*values["H"])

def first(x,values):
    return values["c0"]*np.cosh(np.sqrt((values["k"]/values["D"]))*x)/np.cosh(np.sqrt((values["k"]/values["D"]))*values["H"])

input_file = "input/input.txt"
values = read_input_file(input_file)

isolated=f"output/reaction/isolated/{values['D']:.6f}_{values['k']:.6f}_{values['a']:.6f}_{values['H']:.6f}_{values['c0']:.6f}.txt"

subprocess.run(["run_reactor_isolated.exe"], check=True)
x,c,dc=np.loadtxt(isolated,unpack=True)

plt.figure(figsize=(8, 5))
plt.plot(x, c, linestyle="-", color="royalblue", label="numerical")
plt.plot(x, first(x,values), linestyle="--", color="darkorange", label="analytical")
plt.xlabel("x")
plt.ylabel("c")
plt.title("Priemrjava")
plt.legend()
plt.show()

