#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 09:33:58 2018

@author: janv4
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Legendre
import numpy as np



Ni=200
xi=np.linspace(-1,1,Ni)
yi=np.zeros((Ni))

l=2
for j in range(0,l+1):
    for i in range(0,Ni):
        yi[i] = Legendre.AssociatedLegendre(l,j,xi[i])
    plt.plot(xi,yi,label='P$_{'+str(l)+'}^{'+str(j)+'}$')    


#plt.plot(xi,AssociatedLegendre(2,2,xi),label='P2')
#plt.plot(xi,AssociatedLegendre(3,2,xi),label='P3')
#plt.plot(xi,AssociatedLegendre(4,2,xi),label='P4')
#plt.plot(xi,AssociatedLegendre(5,2,xi),label='P5')
plt.xlabel('x')
plt.ylabel('$P_n(x)$')
plt.legend()
plt.savefig('LegendrePolyPlots.png')
plt.show()

