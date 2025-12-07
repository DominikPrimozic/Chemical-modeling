# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 19:40:53 2025

@author: domin
"""
import numpy as np
import matplotlib.pyplot as plt
import os

flux=f"output/reaction/flux/0.000001_1.000000_1.000000_300.000000_100.000000.txt"
isolated=f"output/reaction/isolated/0.000001_0.001000_1.000000_0.100000_1.000000.txt"
profile=f"output/reaction/profile/0.000001_1.000000_1.000000_300.000000_100.000000.txt"

x,c,dc=np.loadtxt(isolated,unpack=True)

plt.plot(x,c)
plt.show()

x,c,dc,T,dT=np.loadtxt(flux,unpack=True)

plt.plot(x,c)
#plt.plot(x,T)
plt.show()

x,c,dc=np.loadtxt(profile,unpack=True)

plt.plot(x,c)
plt.show()