#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a module of code that can assist in computing the Legendre 
expansion coefficients.
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import Legendre
import sys



#============================
def Ef(muc,A,Ei):
    alpha = ((A-1)/(A+1))**2.0
    
    return Ei*( (1+alpha) + (1-alpha)*muc )/2

#============================
def Ei(mu,A,Ef):
    alpha = ((A-1)/(A+1))**2.0
    
    return Ef*2/( (1+alpha) + (1-alpha)*mu )

#============================ Probability scattering Energy
def PEftoEi(A,Ei,Ef):
    alpha = ((A-1)/(A+1))**2.0
    
    if (  (Ef<(alpha*Ei)) or (Ef>Ei)  ):
        return 1.0e-10
    
    return (1/(1-alpha)/Ei)

#============================
def ThetaL(A,thetac):
    if (((1/A) + math.cos(thetac))<0):
        return (math.atan(math.sin(thetac)/((1/A) + math.cos(thetac))))+\
               math.pi
    return (math.atan(math.sin(thetac)/((1/A) + math.cos(thetac))))

#============================
def ThetaC(thetaL,A):
    B = 1/A
    if (thetaL == (math.pi/2)):
        C=1
    else:
        C = min(math.tan(thetaL),1.0e6)
    
    root = (1-B**2+1/(C**2))
    
    x1 = (-2*B*C**2 + 2*C*math.sqrt(root))/2/(C**2+1)
    #x2 = (-2*B*C**2 - 2*C*math.sqrt(root))/2/(C**2+1)
    
    #Safety catches
    if (x1>1.0):
        return math.acos(1)
    if (x1<-1.0):
        return math.acos(-1)
    
    return math.acos(x1)

#============================ Probability scattering mu
def Pmu(mu,A):
    thetaLA = math.acos(mu-0.0000001)
    thetaLB = math.acos(mu+0.0000001)
    thetacA = ThetaC(thetaLA,A)
    thetacB = ThetaC(thetaLB,A)
    CPLA    = 0.5-0.5*math.cos(thetacA)
    CPLB    = 0.5-0.5*math.cos(thetacB)
    return -(CPLB - CPLA)/0.0000002

#============================ Kernel
def Kernel(mu,gprime,g,Eg,A,Ng=300):
    #=== Bin boundaries
    Eiupp = Eg[gprime]
    Eilow = Eg[gprime+1]
    Efupp = Eg[g]
    Eflow = Eg[g+1]
    
    dEi = (Eiupp - Eilow)/Ng
    binWidth = (Efupp-Eflow)
    
    sumprobs=0
    thetaL = math.acos(mu)
    thetac = ThetaC(thetaL,A)
    muc = math.cos(thetac)
    for iE in range(0,Ng):
        Ein = Eilow + dEi/2 + dEi*iE        
        Eout=Ef(muc,A,Ein)
        
        if ((Eout<=Efupp) and (Eout>=Eflow)):
            sumprobs = sumprobs + 1/Ng
    
    return sumprobs*Pmu(mu,A)

#============================ Find KernelL(E'toE)
def KernelL(ell,gprime,g,Eg,A):
    groupprob=0
    Np=800
    mu=np.linspace(-0.9999,0.9999,Np)
    dmu = (0.9999*2)/Np
    print("Integrating group %d to %d moment %d" %(gprime,g,ell))
    for i in range(0,Np):
        groupprob = groupprob + Kernel(mu[i],gprime,g,Eg,A)* \
             Legendre.Legendre(ell,mu[i])*dmu
    
    return groupprob

#============================ Expansion function
def expansion(L,mu,KL):
    value=0
    for ell in range(0,L+1):
        value = value +((2*ell+1)/2)*KL[ell]*Legendre.Legendre(ell,mu)
    return value
