# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 20:44:14 2025

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 22:23:07 2025

@author: domin
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

folder_path = 'output/porous/concentrations'

# Get a sorted list of filenames
filenames = sorted(os.listdir(folder_path), key=lambda x: int(x.replace('.txt', '')))

i = 0
for filename in filenames:
    if filename.endswith('.txt'):
        if i % 2 == 0:
            print(i)
            concentration = int(filename.replace('.txt', ''))
            file_path = os.path.join(folder_path, filename)
            _, _, c = np.loadtxt(file_path, unpack=True)

            plt.plot(c, label=concentration)
        i += 1

plt.xlabel("Časovni koraki")
plt.ylabel("Koncentracija [delci/ploščina]")
plt.legend(title="Delci medija", loc='upper left', bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.show()

plt.savefig("sprememba_cc_s_spremembo_števila_internih.png")
