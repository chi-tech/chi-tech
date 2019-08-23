import numpy as np

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt

cmapjet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0.0,vmax=1.0)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmapjet)

import ZTET

default_points = []
default_points.append(np.array([0.0,0.0,0.0]))
default_points.append(np.array([0.0,1.0,0.0]))
default_points.append(np.array([1.0,1.0,0.0]))
default_points.append(np.array([1.0,0.0,0.0]))
default_points.append(np.array([0.0,0.0,1.0]))
default_points.append(np.array([0.0,1.0,1.0]))
default_points.append(np.array([1.0,1.0,1.0]))
default_points.append(np.array([1.0,0.0,1.0]))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class QuadFace:
  def __init__(self,p0,p1,p2,p3):
    self.v_indices = np.array([p0,p1,p2,p3],dtype=int)
    self.edges = []
    self.edges.append(np.array([p0,p1],dtype=int))
    self.edges.append(np.array([p1,p2],dtype=int))
    self.edges.append(np.array([p2,p3],dtype=int))
    self.edges.append(np.array([p3,p0],dtype=int))

    self.fc=np.zeros(3)
    self.sides = []
    self.sides_indices = []

  def ComputeFaceCenter(self,nodes):
    for p in range(0,4):
      self.fc += nodes[self.v_indices[p]]
    self.fc /= 4

  def ComputeSides(self,nodes,vc):
    self.ComputeFaceCenter(nodes)
    for e in range(0,len(self.edges)):
      p0_i = self.edges[e][0]
      p1_i = self.edges[e][1]

      p0 = nodes[p0_i]
      p1 = nodes[p1_i]
      p2 = vc
      p3 = self.fc

      self.sides.append(ZTET.Tet(p0,p1,p2,p3))
      self.sides_indices.append([p0_i,p1_i])




  def PlotFace(self,ax,nodes):
    for e in range(0,len(self.edges)):
      points = np.zeros((2,3))
      points[0,:] = nodes[self.edges[e][0]]
      points[1,:] = nodes[self.edges[e][1]]

      ax.plot(points[:,0],points[:,1],points[:,2],'k')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class HexaHedral:
  def __init__(self,points=default_points):
    self.nodes = points
    self.faces = []
    self.faces.append(QuadFace(0,4,5,1))
    self.faces.append(QuadFace(2,1,5,6))
    self.faces.append(QuadFace(2,6,7,3))
    self.faces.append(QuadFace(0,3,7,4))
    self.faces.append(QuadFace(0,1,2,3))
    self.faces.append(QuadFace(7,6,5,4))

    # Compute cell center
    self.vc = np.zeros(3)
    for p in range(0,len(points)):
      self.vc += points[p]
    self.vc /= len(points)

    # Compute sides
    for f in range(0,len(self.faces)):
      self.faces[f].ComputeSides(self.nodes,self.vc)

    # Compute node_side map
    self.node_maps = []
    for n in range(0,len(self.nodes)):
      FEnodeMap = []
      for f in range(0,len(self.faces)):
        FEnodeFaceMap = []
        for s in range(0,len(self.faces[f].sides)):
          side_indices = self.faces[f].sides_indices[s]

          if   (side_indices[0] == n):
            sideMap = 0
          elif (side_indices[1] == n):
            sideMap = 1
          else:
            sideMap = -1

          FEnodeFaceMap.append(sideMap)
        FEnodeMap.append(FEnodeFaceMap)
      self.node_maps.append(FEnodeMap)




  def PlotOutline(self,ax):
    for f in range(0,len(self.faces)):
      self.faces[f].PlotFace(ax,self.nodes)

  def PlotSides(self,ax):
    for f in range(0,len(self.faces)):
      for s in range(0,len(self.faces[f].sides)):
        self.faces[f].sides[s].PlotOutline(ax,'g')

  def Varphi(self,i,x,y,z):
    for f in range(0,len(self.faces)):
      for s in range(0,len(self.faces[f].sides)):
        if (self.faces[f].sides[s].IsPointInside(x,y,z)):
          varphi = 0.0
          varphi_i  = 0.0
          varphi_fc = 0.0
          varphi_cc = 0.0

          if (self.node_maps[i][f][s] >= 0):
            mapping=self.node_maps[i][f][s]
            varphi_i=self.faces[f].sides[s].Varphi(mapping,x,y,z)

          for v in range(0,len(self.faces[f].v_indices)):
            if (self.faces[f].v_indices[v]== i) :
              varphi_fc=self.faces[f].sides[s].Varphi(3,x,y,z)


          varphi_cc=self.faces[f].sides[s].Varphi(2,x,y,z)


          varphi=varphi_i+(1/4)*varphi_fc+(1/8)*varphi_cc
          return varphi

  def GradVarphi(self,i,x,y,z):
    grad_varphi=np.zeros(3)
    for f in range(0,len(self.faces)):
      for s in range(0,len(self.faces[f].sides)):
        if (self.faces[f].sides[s].IsPointInside(x,y,z)):
          grad_varphi    = np.zeros(3)
          grad_varphi_i  = np.zeros(3)
          grad_varphi_fc = np.zeros(3)
          grad_varphi_cc = np.zeros(3)

          if (self.node_maps[i][f][s] >= 0):
            mapping=self.node_maps[i][f][s]
            grad_varphi_i=self.faces[f].sides[s].GradVarphi(mapping,x,y,z)

          for v in range(0,len(self.faces[f].v_indices)):
            if (self.faces[f].v_indices[v]== i) :
              grad_varphi_fc=self.faces[f].sides[s].GradVarphi(3,x,y,z)

          grad_varphi_cc=self.faces[f].sides[s].GradVarphi(2,x,y,z)

          grad_varphi=grad_varphi_i+(1/4)*grad_varphi_fc+(1/8)*grad_varphi_cc
          return grad_varphi
    return grad_varphi

  def NumIntV_Varphi(self,i,xr,yr,zr):
    Np=len(xr)
    sum = 0.0
    for f in range(0,len(self.faces)):
      for s in range(0,len(self.faces[f].sides)):
        xyz=np.zeros((Np,3))
        col=[]

        for p in range(0,Np):
          xyz[p,:]=self.faces[f].sides[s].GetXYZ(xr[p],yr[p],zr[p])
          mapping=self.node_maps[i][f][s]
          varphi_i=0.0
          varphi_fc=0.0
          varphi_cc=0.0

          if (self.node_maps[i][f][s]>=0):
            varphi_i=self.faces[f].sides[s].Varphi(mapping,
                                                   xyz[p,0],
                                                   xyz[p,1],
                                                   xyz[p,2])
          for v in range(0,len(self.faces[f].v_indices)):
            if (self.faces[f].v_indices[v]==i):
              varphi_fc=self.faces[f].sides[s].Varphi(2,
                                                      xyz[p,0],
                                                      xyz[p,1],
                                                      xyz[p,2])
          varphi_cc=self.faces[f].sides[s].Varphi(3,
                                                  xyz[p,0],
                                                  xyz[p,1],
                                                  xyz[p,2])

          varphi=varphi_i+(1/4)*varphi_fc+(1/8)*varphi_cc

          sum += varphi*self.faces[f].sides[s].detJ
          # sum += varphi

    return sum

  def ScatterVarphi(self,i,xr,yr,zr,ax):
    Np = len(xr)
    for f in range(0,len(self.faces)):
      for s in range(0,len(self.faces[f].sides)):
        xyz=np.zeros((Np,3))
        col=[]

        for p in range(0,Np):
          xyz[p,:]=self.faces[f].sides[s].GetXYZ(xr[p],yr[p],zr[p])
          mapping=self.node_maps[i][f][s]
          varphi_i = 0.0
          varphi_fc= 0.0
          varphi_cc= 0.0

          if (self.node_maps[i][f][s]>=0):
            varphi_i=self.faces[f].sides[s].Varphi(mapping,
                                                   xyz[p,0],
                                                   xyz[p,1],
                                                   xyz[p,2])
          for v in range(0,len(self.faces[f].v_indices)):
            if (self.faces[f].v_indices[v]== i):
              varphi_fc=self.faces[f].sides[s].Varphi(2,
                                                      xyz[p,0],
                                                      xyz[p,1],
                                                      xyz[p,2])
          varphi_cc=self.faces[f].sides[s].Varphi(3,
                                                  xyz[p,0],
                                                  xyz[p,1],
                                                  xyz[p,2])

          varphi=varphi_i+(1/4)*varphi_fc+(1/8)*varphi_cc

          col.append(scalarMap.to_rgba(varphi))

        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c=col,cmap='jet',depthshade=False)





