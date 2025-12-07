# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 18:30:15 2025

@author: domin
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

path="output/concnetration.txt"

a,b,c=np.loadtxt(path, unpack=True)

plt.plot(a, label = "A")
plt.plot(b, label = "B")
plt.plot(c, label = "C")

plt.xlabel("Časovni koraki")
plt.ylabel("Koncentracija [delci/ploščina]")
plt.legend()

plt.savefig("40_0.25_40_0.25_long.png")