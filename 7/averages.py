# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 16:00:01 2025

@author: domin
"""

import numpy as np
import matplotlib.pyplot as plt

path = "output/MDLJ/4/thermo.txt"

t, T, K, U, E, p = np.loadtxt(path, skiprows=1, unpack=True)

# Temperature vs Time
plt.plot(t, T)
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.axhline(y=0.75, color='g', linestyle='--')
plt.savefig('100_0.4_T.png')  # Save the plot as a PNG file
plt.show()
plt.close()  # Close the current plot

# Potential Energy vs Time
plt.plot(t, U)
plt.xlabel('Time')
plt.ylabel('Potential Energy (U)')
# Add dotted orange line for average
avg_U = np.mean(U)
plt.axhline(y=avg_U, color='orange', linestyle=':')
plt.savefig('100_0.4_U.png')  # Save the plot as a PNG file
plt.show()
plt.close()  # Close the current plot

# Kinetic Energy vs Time
plt.plot(t, K)
plt.xlabel('Time')
plt.ylabel('Kinetic Energy (K)')
# Add dotted orange line for average
avg_K = np.mean(K)
plt.axhline(y=avg_K, color='orange', linestyle=':')
plt.savefig('100_0.4_K.png')  # Save the plot as a PNG file
plt.show()
plt.close()  # Close the current plot

# Total Energy vs Time
plt.plot(t, E)
plt.xlabel('Time')
plt.ylabel('Total Energy (E)')
# Add dotted orange line for average
avg_E = np.mean(E)
plt.axhline(y=avg_E, color='orange', linestyle=':')
plt.savefig('100_0.4_E.png')  # Save the plot as a PNG file
plt.show()
plt.close()  # Close the current plot

# Pressure vs Time
plt.plot(t, p)
plt.xlabel('Time')
plt.ylabel('Pressure (p)')
# Add dotted orange line for average
avg_p = np.mean(p)
plt.axhline(y=avg_p, color='orange', linestyle=':')
plt.savefig('100_0.4_P.png')  # Save the plot as a PNG file
plt.show()
plt.close()  # Close the current plot

