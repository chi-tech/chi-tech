import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Legendre
import numpy as np



Ni=200
xi=np.linspace(-1,1,Ni)
yi=np.zeros((Ni))


plt.plot([-1,1],[Legendre.Legendre(0,-1),Legendre.Legendre(0,1)],label='$P_0(x)$')
plt.plot(xi,Legendre.Legendre(1,xi),label='$P_1(x)$')
plt.plot(xi,Legendre.Legendre(2,xi),label='$P_2(x)$')
plt.plot(xi,Legendre.Legendre(3,xi),label='$P_3(x)$')
plt.plot(xi,Legendre.Legendre(4,xi),label='$P_4(x)$')
plt.plot(xi,Legendre.Legendre(5,xi),label='$P_5(x)$')
plt.xlabel('x')
plt.ylabel('$P_n(x)$')
plt.legend()
plt.savefig('LegendrePolyPlots.png')
plt.show()

