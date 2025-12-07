# -*- coding: utf-8 -*-
"""
Created on Thu May 22 09:14:22 2025

@author: domin
"""

import imageio.v2 as imageio
import glob

frames = []
for filename in sorted(glob.glob("CA/frame_*.ppm")):
    frames.append(imageio.imread(filename))

imageio.mimsave("velika.gif", frames, duration=0.15, loop=0)
