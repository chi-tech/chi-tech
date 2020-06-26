#!/anaconda3/bin/python3
import numpy as np 
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as tck

from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

import mpl_toolkits.mplot3d.art3d as art3d

# ######################################################### Class def
class Quad:
    def __init__(self, vertices,weight=1.0,loud=False):
        #Make edges
        NV = len(vertices)
        edges = []
        edgecentres = []
        for i in range(0,NV):
            if i<(NV-1):
                edges.append([vertices[i],vertices[i+1]])
                edgecentres.append(0.5*(vertices[i]+vertices[i+1]))
            else:
                edges.append([vertices[i],vertices[0]])
                edgecentres.append(0.5*(vertices[i]+vertices[0]))

        #Compute 2d centroid
        centroid2d = np.zeros(2)
        for v in vertices:
            centroid2d += v 
        centroid2d /= len(vertices)

        if loud:
            print(centroid2d)
            print(vertices)

        #Compute 3d vertices
        vertices3d = []
        for v in vertices:
            varphi = v[0]
            theta  = v[1]
            vertices3d.append([math.sin(theta)*math.cos(varphi),
                               math.sin(theta)*math.sin(varphi),
                               math.cos(theta)])

        #Compute 3d centroid
        centroid3d = np.zeros(3)
        for v in vertices3d:
            centroid3d += v
        centroid3d /= len(vertices3d)

        self.edges2d = edges
        self.edgecenters2d = edgecentres
        self.vertices2d = vertices
        self.centroid2d = centroid2d
        self.vertices3d = vertices3d
        self.centroid3d = centroid3d
        self.weight = weight

# ######################################################### 
def SplitElement(cell):
    new_cells = []
    for k in range(0,4):
        ref_back_edge = 3 if (k==0) else k-1
        new_verts = [cell.vertices2d[k],
                     cell.edgecenters2d[k],
                     cell.centroid2d,
                     cell.edgecenters2d[ref_back_edge]]

        new_cell = Quad(new_verts,weight=cell.weight/4.0)             
        new_cells.append(new_cell)

    return new_cells

# #########################################################
def RefineInCone(in_elements,varphi,theta,cone_radius):
    num_elements = len(in_elements)

    new_element_set = []

    num_cells_refined = 0

    ref_vec = np.array([math.sin(theta)*math.cos(varphi),
                        math.sin(theta)*math.sin(varphi),
                        math.cos(theta)])

    for k in range(0,num_elements):
        cell = in_elements[k]

        if (np.dot(ref_vec,cell.centroid3d) > math.cos(cone_radius)):
            new_cells = SplitElement(cell)
            for new_cell in new_cells:
                new_element_set.append(new_cell)
            num_cells_refined += 1
        else:
            new_element_set.append(cell)

    print("Total number of angles after refinement: ",len(new_element_set))
    return new_element_set

# #########################################################
def RefineInConeAboutVec(in_elements,ref_vec,cone_radius):
    num_elements = len(in_elements)

    new_element_set = []

    num_cells_refined = 0

    for k in range(0,num_elements):
        cell = in_elements[k]

        if (np.dot(ref_vec,cell.centroid3d) > math.cos(cone_radius)):
            new_cells = SplitElement(cell)
            for new_cell in new_cells:
                new_element_set.append(new_cell)
            num_cells_refined += 1
        else:
            new_element_set.append(cell)

    print("Total number of angles after refinement: ",len(new_element_set))
    return new_element_set

# ========================================== Define elements
# kNa = 3
# kNp = 3

# dvarphi = 2.0*math.pi/(4*kNa)
# dtheta  = math.pi/(2*kNp)

# elements = []
# for i in range(0,4*kNa):
#     for j in range(0,2*kNp):
#         # Define vertices
#         vertices2d = []
#         vertices2d.append(np.array([(i  )*dvarphi,(j  )*dtheta]))
#         vertices2d.append(np.array([(i+1)*dvarphi,(j  )*dtheta]))
#         vertices2d.append(np.array([(i+1)*dvarphi,(j+1)*dtheta]))
#         vertices2d.append(np.array([(i  )*dvarphi,(j+1)*dtheta]))

#         elements.append(Quad(vertices2d))

azi_angles = [0.261799,
            0.785398,
            1.309,
            1.8326,
            2.35619,
            2.87979,
            3.40339,
            3.92699,
            4.45059,
            4.97419,
            5.49779,
            6.02139]
pol_angles = [0.369607,
            0.848367,
            1.32985,
            1.81174,
            2.29323,
            2.77199]
weights = [8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02,
    8.971E-02,
    1.889E-01,
    2.450E-01,
    2.450E-01,
    1.889E-01,
    8.971E-02]

kNa = len(azi_angles)
kNp = len(pol_angles)

azic = []
azic.append(0.0)
for i in range(0,kNa):
    if (i<(kNa-1)):
        azic.append(0.5*(azi_angles[i]+azi_angles[i+1]))
    else:
        azic.append(2.0*math.pi)

polc = []
polc.append(0.0)
for i in range(0,kNp):
    if (i<(kNp-1)):
        polc.append(0.5*(pol_angles[i]+pol_angles[i+1]))
    else:
        polc.append(math.pi)

elements = []
k=0
for i in range(0,kNa):
    for j in range(0,kNp):
        # Define vertices
        vertices2d = []
        vertices2d.append(np.array([azic[i  ]  , polc[j  ] ]))
        vertices2d.append(np.array([azic[i+1]  , polc[j  ] ]))
        vertices2d.append(np.array([azic[i+1]  , polc[j+1] ]))
        vertices2d.append(np.array([azic[i  ]  , polc[j+1] ]))

        elements.append(Quad(vertices2d,weight=weights[k]))
        k+=1


roi_center = np.array([65.5,85.5,0.5])
src_center = np.array([10.5,10.5,79.5])
ref_vector = roi_center - src_center
ref_vector = ref_vector/np.linalg.norm(ref_vector)

elements = RefineInConeAboutVec(elements,ref_vector,3*math.pi/12)
elements = RefineInConeAboutVec(elements,ref_vector,2*math.pi/12)
elements = RefineInConeAboutVec(elements,ref_vector,2*math.pi/12)
elements = RefineInConeAboutVec(elements,ref_vector,1*math.pi/12)

# =============================================== Printing custom quadratures
print("azi_angles={",end='')
i = 0
j = 0
for cell in elements:
    i += 1
    j += 1

    print("{:.8f}".format(cell.centroid2d[0]),end='')
    if j != len(elements):
        print(", ",end='')

    if i == 8:
        i=0
        print("\n            ",end='')
print("}")

print("pol_angles={",end='')
i = 0
j = 0
for cell in elements:
    i += 1
    j += 1

    print("{:.8f}".format(cell.centroid2d[1]),end='')
    if j != len(elements):
        print(", ",end='')

    if i == 8:
        i=0
        print("\n            ",end='')
print("}")

total_weight = 0.0
print("weights   ={",end='')
i = 0
j = 0
for cell in elements:
    i += 1
    j += 1

    total_weight += cell.weight

    print("{:.8f}".format(cell.weight),end='')
    if j != len(elements):
        print(", ",end='')

    if i == 8:
        i=0
        print("\n            ",end='')
print("}")
print(total_weight,total_weight/math.pi/4.0)

# =============================================== Create plotting entities
points2d = [[],[]]
points3d = [[],[],[]]
elements_2dcolls = []
elements_3dcolls = []
for cell in elements:
    vertices2d = []
    for v in cell.vertices2d:
        vertices2d.append(v/math.pi)
    poly = Polygon(vertices2d,closed=True,linestyle='-',fill=False)
    patch = []
    patch.append(poly)
    coll = PatchCollection(patch)
    coll.set_facecolor([0,0,0,0])
    coll.set_edgecolor([0,0,0,1])
    coll.set_linewidth(0.3)

    elements_2dcolls.append(coll)

    polyh = art3d.Poly3DCollection([cell.vertices3d])
    # w = cell.weight/max(weights)
    w = 1.0
    polyh.set_facecolor((w,w,w,1.0))
    polyh.set_edgecolor((0.0,0.0,0.0,1.0))
    polyh.set_linewidth(0.3)
    polyh.set_alpha(1.0)
    polyh.set_clip_on(True)
    elements_3dcolls.append(polyh)


    points2d[0].append(cell.centroid2d[0]/math.pi)
    points2d[1].append(cell.centroid2d[1]/math.pi)

    points3d[0].append(cell.centroid3d[0]*1.0)
    points3d[1].append(cell.centroid3d[1]*1.0)
    points3d[2].append(cell.centroid3d[2]*1.0)

# =============================================== Plotting
fig = plt.figure(figsize=(18,8))
ax2d = fig.add_subplot(121)
ax2d.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax2d.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax2d.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax2d.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax2d.axes.set_xlim([0.0,2.0])
ax2d.axes.set_ylim([0.0,1.0])

ax2d.scatter(points2d[:][0], points2d[:][1],1.0,c='k')
for coll in elements_2dcolls:
    ax2d.add_collection(coll)


ax3d = fig.add_subplot(122, projection='3d')


ax3d.scatter(points3d[:][0],points3d[:][1],points3d[:][2],s=1.0,depthshade=True,c='w')
for polyh in elements_3dcolls:
    ax3d.add_collection3d(polyh)

plt.show()
