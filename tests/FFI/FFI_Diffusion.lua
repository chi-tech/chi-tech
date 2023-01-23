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
chiVolumeMesherSetMatIDToAll(0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-0.5,0.5)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

D = {1.0,0.01}
Q = {1.0,10.0}
XSa = {1.0,10.0}
function D_coef(i,x,y,z)
    return D[i+1]
end
function Q_ext(i,x,y,z)
    return Q[i+1]
end
function Sigma_a(i,x,y,z)
    return XSa[i+1]
end

-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
eps = 1.0e-6
e_vol = chiLogicalVolumeCreate(RPP,L/2-eps,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-L/2+eps,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,L/2-eps,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-L/2+eps,-1000,1000)
t_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,L/2-eps,1000)
b_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,-L/2+eps)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3
t_bndry = 4
b_bndry = 5

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,t_vol,t_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,b_vol,b_bndry)

--############################################### Add material properties
print(solver)
if (solver == nil) then solver = "CFEM" end

DiffusionSolverCreate = chiCFEMDiffusionSolverCreate
DiffusionSetBCProperty = chiCFEMDiffusionSetBCProperty
if (solver == "DFEM") then
    DiffusionSolverCreate = chiDFEMDiffusionSolverCreate
    DiffusionSetBCProperty = chiDFEMDiffusionSetBCProperty
elseif (solver == "FV") then
    DiffusionSolverCreate = chiFVDiffusionSolverCreate
    DiffusionSetBCProperty = chiFVDiffusionSetBCProperty
end


phys1 = DiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)

DiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"dirichlet",0.0)
DiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"dirichlet",0.0)
DiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"dirichlet",0.0)
DiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"dirichlet",0.0)
DiffusionSetBCProperty(phys1,"boundary_type",t_bndry,"dirichlet",0.0)
DiffusionSetBCProperty(phys1,"boundary_type",b_bndry,"dirichlet",0.0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiSolverGetFieldFunctionList(phys1)
chiExportMultiFieldFunctionToVTK({fflist[1]},solver.."Diff2D_Dirichlet")

--############################################### Slice plot
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

chiFFInterpolationExportPython(slice2)

local handle = io.popen("python ZPFFI00.py")