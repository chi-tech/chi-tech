#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 22:56:08 2018

@author: janv4
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

#============================ Target nucleus and scattering order
A=12
#============================ Define groups
G = 10
Eg = np.linspace(1,0,G+1)
print("Energy group boundaries")
print(Eg)

#============================ Define group 0 xs
Sigs = 15.0


#------------------------------------------------
#============================ Caculate expansion coefficients
L=10
KL=np.zeros((L+1))
gprime = 0
g=1
sumgroups=0
for ell in range(0,L+1):
    KL[ell] = KernelL(ell,gprime,g,Eg,A)
    sys.stdout.flush()



#============================ Discrete form        
Np=100
mu=np.linspace(-0.9999,0.9999,Np)
yexp=np.zeros((Np))
ydis=np.zeros((Np))
sumofdis=0
sumofexp=0
sumovergroupsdis=0
for i in range(0,Np):
    yexp[i] = expansion(L,mu[i],KL)
    ydis[i] = Kernel(mu[i],gprime,g,Eg,A)
    sumofexp=sumofexp+yexp[i]*(2/Np)
    sumofdis=sumofdis+ydis[i]*(2/Np)
    
    for gdes in range(0,G):
        sumovergroupsdis=sumovergroupsdis+ \
                         Kernel(mu[i],gprime,gdes,Eg,A)*(2/Np)
    
print("Sum of discrete=%f" %sumofdis)
print("Sum of expansion=%f" %sumofexp)
print("Sum over all groups (discrete)=%f" %sumovergroupsdis)

plt.figure(1)
#plt.clf()
#plt.plot(mu,ydis,label='Discrete')
plt.plot(mu,yexp,label='Expansion L=' + str(L))
plt.xlim(-1,1)
plt.xlabel(r'$\mu=\cos\theta_L$')
plt.ylabel(r'$K(\mu,E_'+str(gprime)+'->E_'+str(g)+')$')
plt.legend()

plt.show()
#------------------------------------------------



'''
L=0
gprime = 0
sumgroups=0
Np=800
mu=np.linspace(-0.9999,0.9999,Np)
dmu = (0.9999*2)/Np
for g in range(0,G):
    groupprob=0
    print("Integrating group %d to %d" %(gprime,g),end='')
    for i in range(0,Np):
        groupprob = groupprob + Kernel(mu[i],gprime,g,Eg,A)* \
             Legendre.Legendre(L,mu[i])*dmu
    sumgroups = sumgroups + groupprob
    print(" Group Prob=%f" %groupprob)
  
print("Sum of probs=%f" %sumgroups)
'''


'''
#============================ Integrate Discrete Kernel over angle
Np = 100*1
Na = 200*1
polar = np.linspace(0.0001,math.pi*0.9999,Np)
azimu = np.linspace(0.0001,2*math.pi*0.99999,Na)
dtheta = (math.pi)/Np
dvarphi= (math.pi*2)/Na

nref = np.array([1.0,0.0,0.0])
ndir = np.array([0.0,0.0,0.0])

xt  = np.zeros((Np*Na))
yt  = np.zeros((Np*Na))
zt  = np.zeros((Np*Na))
prob       = np.zeros((Np*Na))
maxprob=0.0
sumprob=0.0
k=-1
gprime=0
for i in range(0,Na):
    print(i)
    for j in range(0,Np):
        varphi = azimu[i]
        theta  = polar[j]
        
        ndir[0] = math.sin(theta)*math.cos(varphi)
        ndir[1] = math.sin(theta)*math.sin(varphi)
        ndir[2] = math.cos(theta)
        
        mu = np.dot(ndir,nref)
        
        for gdes in range(0,G):
            sumprob=sumprob+ \
              Kernel(mu,gprime,gdes,Eg,A)*math.sin(theta)*dtheta*dvarphi
    
        k=k+1
        xt[k]  = ndir[0]
        yt[k]  = ndir[1]
        zt[k]  = ndir[2]
        prob[k]       = Pmu(mu,A)
        sumprob = sumprob+math.sin(theta)*dtheta*dvarphi
        
        if (prob[k]>maxprob):
            maxprob=prob[k]
           
            
print("Max prob=%f" %maxprob)
print("Avg prob=%f" %(np.sum(prob)/Np/Na))
print("Sum prob=%f" %(sumprob))     
'''

'''
#============================ Integrate P(mu)
Np=50
mu = np.linspace(-0.9999,0.9999,Np)
dmu = (2*0.9999)/Np
probsum=0

for i in range(0,Np):
    probsum = probsum + Pmu(mu[i],A)*dmu
    
    
print(probsum)
'''

'''
plt.figure(1)
plt.clf()
plt.plot(muL,Eout,label='A='+str(A))
plt.xlim(-1,1)
plt.xlabel(r'$\mu=\cos\theta_L$')
plt.ylabel(r'Eout')
plt.legend()

plt.show()
'''




'''
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
plt.colorbar(sc)
plt.show()
'''
