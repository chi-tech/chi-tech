import numpy as np
import math

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import animation

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ref Tet
RefTet = ZTET.Tet()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arb Tet 1 scaled
tet_points = np.zeros((4,3))

scale = 0.5
tet_points[0,:] = np.array([0.0,0.0,0.0])
tet_points[1,:] = np.array([1.0,0.0,0.])
tet_points[2,:] = np.array([0.0,1.0,0.0])
tet_points[3,:] = np.array([0.0,0.0,1.0])

ArbTet = ZTET.Tet(tet_points[0,:],
                  tet_points[1,:],
                  tet_points[2,:],
                  tet_points[3,:])

p0 = ArbTet.p0*1.0 + ArbTet.cellcenter*0.0
p1 = ArbTet.p0*0.0 + ArbTet.cellcenter*1.0

dof=1
v0 = ArbTet.Varphi(dof,p0[0],p0[1],p0[2])
v1 = ArbTet.Varphi(dof,p1[0],p1[1],p1[2])

dvref= ArbTet.GradVarphi(dof,p1[0],p1[1],p1[2])

dp = p1-p0
dv = v1-v0

print(dv,dp)
print(ArbTet.detJ)


print("Grad: [%g %g %g]" %(dv/dp[0],dv/dp[1],dv/dp[2]))
print("GradR:[%g %g %g]" %(dvref[0],dvref[1],dvref[2]))






# ========================================== Plotting
fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')

ArbTet.PlotOutline(ax)
ArbTet.PlotGrad(0,ax)


plt.show()