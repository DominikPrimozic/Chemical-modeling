# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 18:15:39 2025

@author: domin
"""
import subprocess
import numpy as np
import matplotlib.pyplot as plt


file_path = "input/duffing.txt"
params = {}


with open(file_path, 'r') as file:
    params['a'], params['b'], params['d'], params['g'], params['w'] = map(float, file.readline().split())
    params['x0'], params['dx0'] = map(float, file.readline().split())
    params['steps'] = int(file.readline())
    params['dt'] = float(file.readline())
    params['file_output'] = file.readline().strip()


output_file = f"output/duffing/{params['file_output']}.txt"


subprocess.run(["run_duffing.exe"], check=True)



t, x, dx = np.loadtxt(output_file , unpack=True)

plt.figure(figsize=(10, 5))
plt.plot(t, x, label="x(t)")
plt.plot(t, dx, label="dx/dt", linestyle='dashed')
plt.xlabel("Time t")
plt.ylabel("Amplitude")
plt.legend()
plt.title(f"$\\alpha$={params['a']}, $\\beta$={params['b']}, $\\delta$={params['d']}, $\\gamma$={params['g']}, $\\omega$={params['w']}")

plt.text(1.05, 0.5, f"$x(0)$={params['x0']}\n $x'(0)$={params['dx0']}\n $\\Delta t$={params['dt']}",
         horizontalalignment='left', verticalalignment='center',
         transform=plt.gca().transAxes, fontsize=12)

plt.subplots_adjust(right=0.85)

output_image = f"output/duffing/{params['file_output']}.png"
plt.savefig(output_image, dpi=300) 
plt.show()