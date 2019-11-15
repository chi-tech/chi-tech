import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Legendre
import math
import numpy as np


print('Hello')
Ni=200
xi=np.linspace(-1,1,Ni)
yi=np.zeros((Ni))

Np = 8
Na = 8

print("PolarLegendre")
xnp,wnp = Legendre.LegendreRoots(2*Np)
print("Legendre")
xna,wna = Legendre.LegendreRoots(4*Na)

print("Sum of Polar weights=%f" % np.sum(wnp))
print("Sum of Azimuthal weights=%f" % np.sum(wna))

NpNa=Np*Na*4*2
varphi = np.zeros((NpNa))
theta = np.zeros((NpNa))
weight = np.zeros((NpNa))
x = np.zeros((NpNa))
y = np.zeros((NpNa))
z = np.zeros((NpNa))
col = np.zeros((NpNa))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#============================================
def f0(i_varphi,i_theta):
    return 1.0;    
    
def f1(i_varphi,i_theta):
    return 1.0 + math.cos(i_theta)*math.sin(i_varphi);

def f2(i_varphi,i_theta):
    return 1.0 + math.pow(math.cos(i_theta),2.0)*math.sin(i_varphi/2);


#============================================
weightsum=0;
k=0
for i in range(0,Na*4):
    for j in range(0,Np*2):
        varphi[k] = math.pi*xna[i]+math.pi
        #varphi[k] = math.pi*(2*(i+1)-1)/Na/8
        theta[k] = math.acos(xnp[j])
        weight[k] = wna[i]*wnp[j]*math.pi
        #weight[k] = wnp[j]*wna[i]*2
        #print('Angle pair (%f,%f), weight=%f' % (varphi[k]*180/math.pi,theta[k]*180/math.pi,weight[k]))
        xt=math.sin(theta[k])*math.cos(varphi[k])
        yt=math.sin(theta[k])*math.sin(varphi[k])
        zt=math.cos(theta[k])
        weightsum=weightsum+weight[k]*f0(varphi[k],theta[k])
        
        if ((xt>=0) and (yt>=0) and (zt>=0)):
            x[k]=math.sin(theta[k])*math.cos(varphi[k])
            y[k]=math.sin(theta[k])*math.sin(varphi[k])
            z[k]=math.cos(theta[k])
            col[k]=weight[k]
            #print("%d %.3f %.3f %.1f %.1f" %(k,varphi[k], \
            #    theta[k],varphi[k]*180/math.pi,theta[k]*180/math.pi))
            k=k+1

 
sc = ax.scatter(x,y,z,s=col*1000,c=col,cmap='jet') 

 
plt.xlabel("X")
plt.ylabel("Y")
ax.set_zlabel('Z')  
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
plt.colorbar(sc)
plt.show()
print("Integration = %f vs %f " %(weightsum, 4*math.pi+8/3))

