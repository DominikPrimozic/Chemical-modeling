# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 13:26:42 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 18:15:39 2025

@author: domin
"""
import subprocess
import numpy as np
import matplotlib.pyplot as plt

maxx= []

output_file="output/duffing/duffin_1_1_1_1_0"

def xx(maxx,output_file): 
    _, x, _ = np.loadtxt(output_file , unpack=True)
    maxx.append(max(x))

output_file="output/duffing/duffin_1_1_1_1_0.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_1.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_2.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_3.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_4.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_5.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_6.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_7.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_8.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_9.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_10.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_11.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_12.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_13.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_14.txt"
xx(maxx,output_file)
output_file="output/duffing/duffin_1_1_1_1_15.txt"
xx(maxx,output_file)

plt.plot(maxx)
"""

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
"""