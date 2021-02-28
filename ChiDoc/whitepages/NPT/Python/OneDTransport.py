#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 20:32:23 2018

@author: janv4
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import Legendre


L=1.0       #Domain length [0,L]
Ne=10    #Number of elements
dx = L/Ne

#============================ Shape functions
def shape(j,xi):
    if (j==0):
        return -0.5*xi + 0.5
    else:
        return 0.5*xi+0.5
    
def gradxshape(j,xi):
    if (j==0):
        return -0.5
    else:
        return 0.5
    
#============================ Quadrature points
qpoints=2
xn,wn = Legendre.LegendreRoots(qpoints)
print(xn,wn)

#============================ Jacobian
detJ = L/Ne/2

#============================ Left BC
PsiL = 1.0
PsiR = 0.0

#============================ Source function
def Q(x):
    value = math.cos(2*math.pi*x)
    
    return value

#============================ Compute numerical solution
Psi_iter = []
plt.figure(1)
plt.clf()
for c in range(0,Ne):
    A = np.zeros((2,2))
    b = np.zeros((2))
    x = np.zeros((2))
    x[0] = c*dx
    x[1] = (c+1)*dx
    xc = 0.5*dx + c*dx
    
    for q in range(0,qpoints):
        A[0,0] = A[0,0] + wn[q]*shape(0,xn[q])* \
                           gradxshape(0,xn[q])*detJ/detJ
        A[0,1] = A[0,1] + wn[q]*shape(0,xn[q])* \
                           gradxshape(1,xn[q])*detJ/detJ
        A[1,0] = A[1,0] + wn[q]*shape(1,xn[q])* \
                           gradxshape(0,xn[q])*detJ/detJ
        A[1,1] = A[1,1] + wn[q]*shape(1,xn[q])* \
                           gradxshape(1,xn[q])*detJ/detJ         
        b[0] = b[0] + Q(xc)*wn[q]* \
                                shape(0,xn[q])*detJ
        b[1] = b[1] + Q(xc)*wn[q]* \
                                shape(1,xn[q])*detJ
        
    A[0,0] = A[0,0] + 0.5           
    A[1,1] = A[1,1] 
          
    if (c==0):
        b[0] = b[0] + 0.5*PsiL
    else:
        b[0] = b[0] + 0.5*PsiR
        
    Psi_local = np.linalg.solve(A,b)
    
    PsiR = Psi_local[1]
    Psi_iter.append(Psi_local)
    plt.plot(x,Psi_local)
    
    
#============================ Compute analytical solution
x = np.linspace(0,L,200)
y = np.zeros((200))
for i in range(0,200):
    y[i] = (1/2/math.pi)*math.sin(2*math.pi*x[i])+1 
    
    
plt.plot(x,y)    
    
plt.show()
            
            
            
            
            
            
            
            