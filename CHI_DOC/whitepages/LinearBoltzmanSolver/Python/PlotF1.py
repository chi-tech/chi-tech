import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Legendre
import numpy as np
import math
import scipy.special
import random



#============================
def F0(varphi,theta):
    return 1

def F1(varphi,theta):
    return 1.0+0.1*math.cos(varphi*4)

def F2(varphi,theta):
    if (varphi<(7*math.pi/8)):
        return 0.2;
    if (varphi>(9*math.pi/8)):
        return 0.2;
    
    return 1.2*math.cos(varphi*4)+0.2

def F3(varphi,theta):
    if (math.cos(varphi*8)<0):
        return 0.2;
    
    return 1.2*math.cos(varphi*8)+0.2

def F4(varphi,theta):
    if (varphi<(math.pi/4)):
        return 0.2;
    if (varphi>(math.pi/2)):
        return 0.2;
    
    return max(1.2*math.cos(varphi*6)+0.2,0.2)

#============================
'''
def Flm(ell,m,varphi,theta):
    return F3(varphi,theta)*Min1powerm(m)* \
           (Legendre.Yl_m(ell,-m,varphi,theta))
'''           
def Flm(ell,m,varphi,theta):
    return F1(varphi,theta)*Legendre.Min1powerm(0)* \
           (Legendre.Ylmcoeff(ell,m,varphi,theta))
          
def Ylmlm(ell,m,varphi,theta)   :
    return Legendre.Ylmcoeff(ell,m,varphi,theta)* \
           Legendre.Ylmcoeff(ell,m,varphi,theta)         

GLC = Legendre.Quadrature()
GLC.InitializeWithGLC(16,8)


#=========================== Build flm
L =4
k=-1
flm = np.zeros((L*(L+2)+1),dtype=np.complex_)
for ell in range(0,L+1):
    for m in range(-ell,ell+1):
        k=k+1
        #flm[k] = Legendre.RiemannAngLM(F1lm,ell,m)
        
        #print("l=%f, m=%f, flm=" %(ell,m), end='')
        #print(flm[k])
        flm[k] = Legendre.QuadratureIntegrateLM(Flm,GLC,ell,m)
        print("l=%d, m=%d, flm=" %(ell,m), end='')
        print("%.5f" %(flm[k]))

#============================ Build yi1
Ni1=400
varphi1=np.linspace(0,2*math.pi,Ni1)
yi1=np.zeros((Ni1))

for i in range(0,Ni1):
    yi1[i]= 0;
    #v = 0+0j
    v=0
    k=-1;
    for ell in range(0,L+1):
        for m in range(-ell,ell+1):
            k=k+1
            v=v+(flm[k])* \
            (Legendre.Ylmcoeff(ell,m,varphi1[i],math.pi/2))
            yi1[i] = v
    print("Fstar %g" %v)

#============================ Build yi2
Ni2=200
varphi2=np.linspace(0,2*math.pi,Ni2)
yi2=np.zeros((Ni2))
for i in range(0,Ni2):
    yi2[i] = F1(varphi2[i],math.pi/2)

'''
#============================ Testing orthogonality
Ltest=4
for ell in range(0,Ltest+1):
        for m in range(-ell,ell+1):
            print("l=%f, m=%f, int Ylm*Ylm=" %(ell,m), end='')
            print(Legendre.RiemannAngLM(Ylmlm,ell,m))  
'''

#============================ Testing addition theorem
Ltest=12
theta_a = math.pi*random.random();
phi_a   = 2*math.pi*random.random();
theta_b = math.pi*random.random();
phi_b   = 2*math.pi*random.random();

n_a = [math.sin(theta_a)*math.cos(phi_a),
       math.sin(theta_a)*math.sin(phi_a),
       math.cos(theta_a)]
n_b = [math.sin(theta_b)*math.cos(phi_b),
       math.sin(theta_b)*math.sin(phi_b),
       math.cos(theta_b)]
print(np.dot(n_a,n_b))


for ell in range(0,Ltest+1):
        Pl = Legendre.Legendre(ell,np.dot(n_a,n_b))
        sums = 0
        for m in range(-ell,ell+1):
            sums = sums+(4*math.pi/(2*ell+1))* \
                   Legendre.Ylmcoeff(ell,m,phi_a,theta_a)* \
                   Legendre.Ylmcoeff(ell,m,phi_b,theta_b)
            
        # print("Pl - sums=", end='')
        # print(Pl - sums)

'''
plt.clf()
plt.plot(varphi2*180/math.pi,yi2,label=r'$f(\varphi,\theta)$')
plt.plot(varphi1*180/math.pi,yi1,label='Expansion L='+str(L),
        marker='D',fillstyle='none',color='k',linestyle='',
        markersize=3,markeredgewidth=.5)
plt.xlabel(r'$\varphi$ [degrees]')
plt.ylabel(r'f($ \varphi, \theta$ )')
axes = plt.gca()
#axes.set_xlim([xmin,xmax])
axes.set_ylim([0,1.5])
plt.legend()

'''
plt.clf()
ax=plt.subplot(111,projection='polar')

ax.plot(varphi2,yi2,label=r'$f(\varphi,\theta)$')
ax.plot(varphi1,yi1,label='Expansion L='+str(L),
        marker='D',fillstyle='none',color='k',linestyle='',
        markersize=3,markeredgewidth=.5)


plt.legend(loc='upper right',bbox_to_anchor=(1.2,1.1))
ax.set_rmax(1.2)
#ax.set_rticks([-0.5, 1, 1.5])
ax.set_rticks([0, 1, 1.5])

plt.show()


# Nt = 10
# Nv = 10
# dtheta = math.pi/Nt
# dvaphi = math.pi*2/Nv
#
# ell = 7
# m   = 1
#
# for i in range(0,Nt):
#   for j in range(0,Nv):
#     theta = dtheta/2 + i*dtheta
#     vaphi = dvaphi/2 + j*dvaphi
#     Ylm = Legendre.Ylmcoeff(ell,m,vaphi,theta)
#     #Ylm = Legendre.AssociatedLegendre(ell,abs(m),math.cos(theta))
#     print("Ylm = %g t=%g v=%g" %(Ylm,theta,vaphi))

