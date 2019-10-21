#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 22:56:08 2018

@author: janv4
"""
import matplotlib.pyplot as plt
import numpy as np
import math
from pyne import data
from pyne import ace
#print(data.atomic_mass('U235'))

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
    return -(CPLB - CPLA)/0.0000001

#============================ Kernel
def Kernel(mu,Ei,Ef,A):
    return 0

#============================
# pathtoH = '/Users/janv4/Desktop/Projects/MCNPDATA/MCNP_DATA'+ \
#           '/DVD-3/MORE_DATA_620-B/xdata/endf71x/H/1001.710nc'
#
# libFile = ace.Library(pathtoH)
# libFile.read()
#
# H = libFile.tables['1001.710nc']
#
# elastic = H.reactions[2]
# radcap = H.reactions[102]

#============================ Angle-to-Energy


#============================ Build CPF
Ni=50
thetac = np.linspace(0.0001,math.pi*0.99999,Ni)
Pc     = np.zeros((Ni+1))
thetaL = np.zeros((Ni+1))
CPL     = np.zeros((Ni))
muL    = np.zeros((Ni))
Einit  = 1.0
Efinal = np.zeros((Ni))
A      = 208
Pc[Ni] = 0
thetaL[Ni]=math.pi

for i in range(0,Ni):
    Pc[i]     = (math.cos(thetac[i])+1)/2
    thetaL[i] = ThetaL(A,thetac[i])
    muL[i]    = math.cos(thetaL[i])
    Efinal[i] = Ef(math.cos(thetac[i]),A,Einit)
    CPL[i]     = (1-math.cos(ThetaC(thetac[i],A)))/2

#============================ PDF
Np=1000
mu = np.linspace(-0.99999,0.99999,Np)
PL = np.zeros((Np))
sumPL=0
dmu=2/Np
for i in range(0,Np):
    thetaLA = math.acos(mu[i]-0.0000001)
    thetaLB = math.acos(mu[i]+0.0000001)
    thetacA = ThetaC(thetaLA,A)
    thetacB = ThetaC(thetaLB,A)
    CPLA    = 0.5-0.5*math.cos(thetacA)
    CPLB    = 0.5-0.5*math.cos(thetacB)
    PL[i] = -(CPLB - CPLA)/0.0000001
    sumPL = sumPL +PL[i]*dmu
    
print(sumPL)
print(math.acos(-1))



'''
#============================ Plot muL to Energy
plt.figure(1)
plt.clf()

plt.plot(thetac*180/math.pi,thetaL[0:(Ni)]*180/math.pi)
plt.xlabel(r'$\theta_c$ [degrees]')
plt.ylabel(r'$\theta_L$ [degrees]')
'''

'''
plt.figure(2)
plt.clf()
plt.plot(thetac*180/math.pi,Efinal)
plt.xlabel(r'$\theta_c$ [degrees]')
plt.ylabel(r'$E_L final $ [degrees]')
'''

#============================ Plot CPF
plt.figure(3)
plt.clf()
plt.plot(thetac*180/math.pi,CPL,label='A='+str(A))
#plt.plot(thetaL*180/math.pi,Pc,label='A='+str(A)+'(Numerical)',
#         marker='D',linestyle='',color='k',fillstyle='none')
plt.xlim(0,180)
plt.xlabel(r'$\theta_L$ [degrees]')
plt.ylabel(r'Cumulative Probability')
plt.legend()

#============================ Plot CPD
plt.figure(4)
plt.clf()
plt.plot(mu,PL,label='A='+str(A))
plt.xlim(-1,1)
plt.xlabel(r'$\mu=\cos\theta_L$')
plt.ylabel(r'P($\mu$)')
plt.legend()


'''
print(H.reactions.items())
fig1, ax1 = plt.subplots()
ax1.loglog(H.energy,H.sigma_t,label='Total')
ax1.loglog(H.energy,elastic.sigma,label='Elas')
ax1.loglog(H.energy,radcap.sigma,label='Capture')
ax1.legend()
'''



plt.show()