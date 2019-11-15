import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Legendre
import math

#============================ Test functions
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

#============================ Quadrature
F=F4
quadGLL = Legendre.Quadrature()
quadGLL.InitializeWithGLL(16,2)

print("Gauss-Legendre-Legendre Quadrature Calculated integral:" )
print(Legendre.QuadratureIntegrate(F,quadGLL))

quadGLC = Legendre.Quadrature()
quadGLC.InitializeWithGLC(16,2)

print("Gauss-Legendre-Chebyshev Quadrature Calculated integral:" )
print(Legendre.QuadratureIntegrate(F,quadGLC))

print("Expensive Riemann-sum Calculated integral:" )
print(Legendre.RiemannAng(F))
