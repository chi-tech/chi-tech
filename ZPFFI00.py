import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.colors as colors

class Datapoint:
  def __init__(self,x,y,value):
    self.x = x
    self.y = y
    self.value = value

class CellData:
  def __init__(self):
    self.data_points=[]
    self.xy = []
    self.c = []

import ZPFFI01

class BaseDataClass:
  def __init__(self):
    data_object = []
    self.data_object = data_object
    data_object = ZPFFI01.AddData(data_object)

data = BaseDataClass()
print(len(data.data_object))

N = len(data.data_object)

maxavg = 0.0
maxc = 0.0
xycount = 0
for c in range(0,N):
    for i in range(0,np.size(data.data_object[c].c)):
        xycount += 1
        if data.data_object[c].avg>maxavg:
            maxavg = data.data_object[c].avg
        if data.data_object[c].c[i]>maxc:
            maxc = data.data_object[c].c[i]
            
xmax = -9999999.0
xmin =  9999999.0
ymax = -9999999.0
ymin =  9999999.0
zmax = -9999999.0
zmin =  9999999.0
x = np.zeros(xycount)
y = np.zeros(xycount)
z = np.zeros(xycount)
xycount = -1
for c in range(0,N):
    for i in range(0,np.size(data.data_object[c].c)):
        xycount += 1
        x[xycount] = data.data_object[c].xy[i,0]
        y[xycount] = data.data_object[c].xy[i,1]
        z[xycount] = data.data_object[c].c[i]

        if x[xycount]>xmax: xmax = x[xycount]
        if x[xycount]<xmin: xmin = x[xycount]
        if y[xycount]>ymax: ymax = y[xycount]
        if y[xycount]<ymin: ymin = y[xycount]
        if z[xycount]>zmax: zmax = z[xycount]
        if z[xycount]<zmin: zmin = z[xycount]

print("phi_max=%g phi_min=%g" %(zmax,zmin))

fig,ax = plt.subplots(1)
cmapjet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0,vmax=maxavg)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmapjet)

cntr1 = plt.tricontourf(x,y,z,124,cmap=cmapjet)

for c in range(0,N):
    col = scalarMap.to_rgba(data.data_object[c].avg)
    poly = Polygon(data.data_object[c].xy,
                   closed=True,linestyle='-',fill=False)
    patch = []
    patch.append(poly)
    coll = PatchCollection(patch)
    coll.set_facecolor([0,0,0,0])
    coll.set_edgecolor([0,0,0,1])
    coll.set_linewidth(0.3)

    ax.add_collection(coll)

fig.colorbar(cntr1,ax=ax)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
plt.show()
