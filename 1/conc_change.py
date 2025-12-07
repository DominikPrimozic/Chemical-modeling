# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 22:23:07 2025

@author: domin
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

folder_path = 'output/reaction_move'

i=0
# Loop through all the files in the directory
for filename in os.listdir(folder_path):
    if filename.endswith('.txt'):
        # Extract the density from the filename (assuming it is the number before .txt)
        if i%15 == 0:
            print(i)
            concentration = float(filename.replace('.txt', ''))
    
            file_path = os.path.join(folder_path, filename)
            _,_,c=np.loadtxt(file_path, unpack=True)
    
            plt.plot(c, label = concentration)
        i+=1
        

plt.xlabel("Časovni koraki")
plt.ylabel("Koncentracija [delci/ploščina]")
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

# Display the plot
plt.tight_layout()  # Adjust the layout to avoid clipping
plt.show()

plt.savefig("sprememba_cc_s_spremembo_dolžine_premika.png")