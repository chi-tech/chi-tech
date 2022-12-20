--############################################### Setup mesh
chiMeshHandlerCreate()
 
mesh={}
N=100
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
 
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiVolumeMesherExecute();
 
--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

D = {1.0,5.0}
Q = {0.0,1.0}
XSa = {1.0,1.0}
function D_coef(i,x,y,z)
    return D[i+1] -- + x
end
function Q_ext(i,x,y,z)
    return Q[i+1] -- x*x
end
function Sigma_a(i,x,y,z)
    return XSa[i+1]
end

-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.9999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add material properties
--#### DFEM stuff
phys1 = chiDFEMDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)

chiDFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
chiDFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiDFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")
chiDFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 0.0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"square_robin_refl","Flux_Diff")