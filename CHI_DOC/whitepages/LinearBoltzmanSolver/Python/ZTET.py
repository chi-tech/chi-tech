import numpy as np

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):
  def __init__(self,xs,ys,zs,*args,**kwargs):
    FancyArrowPatch.__init__(self,(0,0),(0,0),*args,**kwargs)
    self._verts3d=xs,ys,zs

  def draw(self,renderer):
    xs3d,ys3d,zs3d=self._verts3d
    xs,ys,zs=proj3d.proj_transform(xs3d,ys3d,zs3d,renderer.M)
    self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
    FancyArrowPatch.draw(self,renderer)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Develop reference points
# These are reference points on a reference
# tetrahedron in natural coordinates they can be
# used to define the xy coordinates on any
# arbitrary tet
def GetTetRefPoints(Resolution=10):
  xr = []
  yr = []
  zr = []

  N = Resolution
  ds  = 1/N

  for i in range(0,N):
    # print(i)
    for j in range(0,N):
      for k in range(0,N):
        xi  = ds/2 + i*ds
        eta = ds/2 + j*ds
        zeta= ds/2 + k*ds

        if (xi + eta + zeta)<1.0:
          xr.append(xi)
          yr.append(eta)
          zr.append(zeta)
        else:
          break

  return xr,yr,zr

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RefTet_Shape(i, xi, eta, zeta):
  if (i==0):
    return 1.0 - xi - eta - zeta
  elif (i==1):
    return xi
  elif (i==2):
    return eta
  elif (i==3):
    return zeta
  else:
    print("Error incorrect use of ref tet")
    exit(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RefTet_GradShape_x(i, xi, eta, zeta):
  if (i==0):
    return -1.0
  elif (i==1):
    return  1.0
  elif (i==2):
    return  0.0
  elif (i==3):
    return  0.0
  else:
    print("Error incorrect use of ref tet")
    exit(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RefTet_GradShape_y(i, xi, eta, zeta):
  if (i==0):
    return -1.0
  elif (i==1):
    return  0.0
  elif (i==2):
    return  1.0
  elif (i==3):
    return  0.0
  else:
    print("Error incorrect use of ref tet")
    exit(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RefTet_GradShape_z(i, xi, eta, zeta):
  if (i==0):
    return -1.0
  elif (i==1):
    return  0.0
  elif (i==2):
    return  0.0
  elif (i==3):
    return  1.0
  else:
    print("Error incorrect use of ref tet")
    exit(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Tet:
  def __init__(self,p0=np.array([0,0,0]),
                    p1=np.array([1,0,0]),
                    p2=np.array([0,1,0]),
                    p3=np.array([0,0,1])):
    self.p0 = p0
    self.p1 = p1
    self.p2 = p2
    self.p3 = p3

    self.cellcenter = (p0+p1+p2+p3)/4.0

    self.v01 = p1-p0
    self.v02 = p2-p0
    self.v03 = p3-p0

    self.J = np.zeros((3,3))
    self.J[:,0] = self.v01
    self.J[:,1] = self.v02
    self.J[:,2] = self.v03

    self.Jinv = np.linalg.inv(self.J)

    self.JT    = np.transpose(self.J)
    self.JTinv = np.linalg.inv(self.JT)

    self.detJ = np.linalg.det(self.J)

  def GetXYZ(self,xi,eta,zeta):
    xyz = np.zeros(3)
    xyz = self.p0 + np.matmul(self.J,np.array([xi,eta,zeta]))

    return xyz

  def MapXYZtoRef(self,x,y,z):
    b = np.array([x,y,z]) - self.p0
    xi_eta_zeta = np.matmul(self.Jinv,b)

    return xi_eta_zeta

  def IsPointInside(self,x,y,z):
    xi_eta_zeta = self.MapXYZtoRef(x,y,z)

    xi  = xi_eta_zeta[0]
    eta = xi_eta_zeta[1]
    zeta= xi_eta_zeta[2]

    if ((xi>=-1.0e-4) and (eta>=-1.0e-12) and (zeta>=-1.0e-12) and ((xi+eta+zeta)<=1.0+1.0e-12)):
      return True
    else:
      return False

  def Varphi(self,i,x,y,z):
    xi_eta_zeta = self.MapXYZtoRef(x,y,z)

    value = RefTet_Shape(i,xi_eta_zeta[0],xi_eta_zeta[1],xi_eta_zeta[2])

    return value

  def GradVarphi(self,i,x,y,z):
    xi_eta_zeta = self.MapXYZtoRef(x,y,z)

    value_x = RefTet_GradShape_x(i,xi_eta_zeta[0],xi_eta_zeta[1],xi_eta_zeta[2])
    value_y = RefTet_GradShape_y(i,xi_eta_zeta[0],xi_eta_zeta[1],xi_eta_zeta[2])
    value_z = RefTet_GradShape_z(i,xi_eta_zeta[0],xi_eta_zeta[1],xi_eta_zeta[2])

    ref_grad_xyz = np.array([value_x,value_y,value_z])
    grad_xyz = np.matmul(self.JTinv,ref_grad_xyz)

    return grad_xyz

  def PlotOutline(self,ax,format='k'):
    vab = np.zeros((2,3))
    vab[0,:] = self.p0
    vab[1,:] = self.p1

    line0, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)

    vab=np.zeros((2,3))
    vab[0,:]=self.p0
    vab[1,:]=self.p2

    line1, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)

    vab=np.zeros((2,3))
    vab[0,:]=self.p0
    vab[1,:]=self.p3

    line2, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)

    vab=np.zeros((2,3))
    vab[0,:]=self.p1
    vab[1,:]=self.p2

    line3, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)

    vab=np.zeros((2,3))
    vab[0,:]=self.p2
    vab[1,:]=self.p3

    line4, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)


    vab=np.zeros((2,3))
    vab[0,:]=self.p1
    vab[1,:]=self.p3

    line5, = ax.plot(vab[0:2,0],vab[0:2,1],vab[0:2,2],format)

    return line0,line1,line2,line3,line4,line5

  def PlotGrad(self,i,ax,format='k',incolor="r"):
    p = self.p0
    grad_xyz = self.GradVarphi(i,p[0],p[1],p[2])

    xyz = np.zeros((2,3))
    xyz[0,:] = p
    xyz[1,:] = p + grad_xyz

    ax.plot(xyz[0:2,0],xyz[0:2,1],xyz[0:2,2],format)

    arr = Arrow3D(xyz[0:2,0],xyz[0:2,1],xyz[0:2,2],mutation_scale=20,
                  lw=3,arrowstyle="-|>", color=incolor)

    ax.add_artist(arr)

    return arr
