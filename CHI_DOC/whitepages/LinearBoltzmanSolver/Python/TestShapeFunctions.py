import numpy as np
import math

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib.colors as colors

cmapjet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0.0,vmax=1.0)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmapjet)

# This program integrates the polyhedral shape functions
import ZTET
import ZHEX




Nref=20
xr,yr,zr = ZTET.GetTetRefPoints(Nref)
print("Number of reference points = %d" %len(xr))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arb Tet
tet_points = np.zeros((4,3))

tet_points[0,:] = np.array([0.0,0.0,0.0])
tet_points[1,:] = np.array([1.0,0.0,0.])
tet_points[2,:] = np.array([0.0,1.0,0.0])
tet_points[3,:] = np.array([0.0,0.0,1.0])

ArbTet = ZTET.Tet(tet_points[0,:],
                  tet_points[1,:],
                  tet_points[2,:],
                  tet_points[3,:])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hexahedral
hex_points = []
d = 2/32
zet=0.2
hex_points.append(np.array([0.0,0.0,0.0]))
hex_points.append(np.array([0.0,d  ,0.0]))
hex_points.append(np.array([d  ,d  ,0.0]))
hex_points.append(np.array([d  ,0.0,0.0]))
hex_points.append(np.array([0.0,0.0,zet]))
hex_points.append(np.array([0.0,d  ,zet]))
hex_points.append(np.array([d  ,d  ,zet]))
hex_points.append(np.array([d  ,0.0,zet]))
hex = ZHEX.HexaHedral()


N=100
dw=1/N
pointA = np.array([0.0,1.0,0.0])
pointB = np.array([1.0,1.0,1.0])
xyz=np.zeros((N,3))
varphi = np.zeros(N)
dvdx_n = np.zeros(N)
dvdy_n = np.zeros(N)
dvdz_n = np.zeros(N)
grad_v = []
col=[]
prevpoint = pointB
for i in range(0,N):
  w = dw/2 + i*dw
  point = w*pointA + (1.0-w)*pointB
  xyz[i,:] = point
  varphi[i] = hex.Varphi(1,point[0],point[1],point[2])
  gradvi    = hex.GradVarphi(1,point[0],point[1],point[2])
  grad_v.append(gradvi)
  dv = 0.0
  dx = 0.0
  dy = 0.0
  dz = 0.0
  if i==0:
    dv = varphi[i] - hex.Varphi(1,pointB[0],pointB[1],pointB[2])
    dx = point[0] - pointB[0]
    dy = point[1] - pointB[1]
    dz = point[2] - pointB[2]
  else:
    dv=varphi[i]-varphi[i-1]
    dx=point[0]-prevpoint[0]
    dy=point[1]-prevpoint[1]
    dz=point[2]-prevpoint[2]

  dvdx_n[i]=dv/dx
  dvdy_n[i]=dv/dy
  dvdz_n[i]=dv/dz

  prevpoint = point
  col.append(scalarMap.to_rgba(varphi[i]))

  print("v=%g gradv=[%g %g %g] gradvs=[%g %g %g]"
        %(varphi[i],dvdx_n[i],dvdy_n[i],dvdz_n[i],
          gradvi[0],gradvi[1],gradvi[2]))




# fig = plt.figure(0)
# ax = fig.add_subplot(111, projection='3d')
#
# hex.PlotSides(ax)
#
# # ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c=col,cmap='jet',depthshade=False)
# hex.ScatterVarphi(0,xr,yr,zr,ax)
#
# plt.show()

plt.figure(1)
plt.plot(varphi)

plt.show()

# sum = 0.0
# for i in range(0,8):
#   value = hex.NumIntV_Varphi(i,xr,yr,zr)*((1.0/Nref)*(1.0/Nref)*(1.0/Nref))
#   sum += value
#   print("Integral %d = %g"%(i,value))
#
# print("Sum=%g" %sum)
# print("Ratio=%g" %(sum/0.00078125))

