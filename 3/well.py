# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 09:01:56 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt

x,psi,dpsi=np.loadtxt("jama.txt",unpack=True)

#psi=psi/sum(psi*psi)**0.5

sorted_indices = np.argsort(x)
x= x[sorted_indices]
psi = psi[sorted_indices]
plt.plot(x,psi)