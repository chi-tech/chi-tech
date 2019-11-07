#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  1 13:47:03 2019

@author: janv4
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np



xyz=np.array([
[-1,-1,-1],
[ 1,-1,-1],
[ 1, 1,-1],
[-1, 1,-1]])

xyz[1][0]=2; xyz[2][0]=3;

print(xyz[:,0])
#print(xyz[0:,[0,1]])


fig = plt.figure(1)
ax = fig.add_subplot(111,projection='3d')

ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'k-',linewidth=0.5)

plt.show()
