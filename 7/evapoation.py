# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 21:56:54 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt

path = "output/cluster/S/thermo.txt"

temperatures=[0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5]
energiesU=[]
energiesK=[]
for i in range(1,11):
    path = f"output/cluster/E/{i}/thermo.txt"
    t, T, K, U, E = np.loadtxt(path, skiprows=1, unpack=True)
    energiesU.append(U[-1])
    energiesK.append(K[-1])


plt.plot(temperatures, energiesK, label="T")
plt.plot(temperatures, energiesU, label="U")
plt.xlabel('Temperature')
plt.ylabel('Energy')
plt.legend()
plt.show()
plt.close()  # Close the current plot


