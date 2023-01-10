-- 2D Diffusion test with Vacuum and Reflecting BCs.
-- SDM: PWLD
-- Test: Max-value=2.50000
num_procs = 4
--Also tests integrating with a lua-function




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=32
L=2.0
xmin = -1.0
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

e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.99999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

--e_bndry = chiRegionAddEmptyBoundary(region1);
--w_bndry = chiRegionAddEmptyBoundary(region1);
--n_bndry = chiRegionAddEmptyBoundary(region1);
--s_bndry = chiRegionAddEmptyBoundary(region1);

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)

prop = chiPhysicsMaterialGetProperty(materials[0],SCALAR_VALUE)
if ((prop.is_empty ~=nil) and (not prop.is_empty)) then
    print("Property table populated, value="..tostring(prop.value))
end

--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLD_MIP")
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-8)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,"boundary_type",e_bndry,"vacuum")
chiDiffusionSetProperty(phys1,"boundary_type",w_bndry,"vacuum")
chiDiffusionSetProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiDiffusionSetProperty(phys1,"boundary_type",s_bndry,"reflecting")

--############################################### Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = chiGetFieldFunctionList(phys1)

--############################################### Slice plot
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

--==========================================
xwing=2.0
function IntegrateMaterialVolume(ff_value,mat_id)
    return xwing
end
ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM_LUA,"IntegrateMaterialVolume")
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
print(chiFFInterpolationGetValue(curffi))
--==========================================

--############################################### Exports
if (master_export == nil) then
    chiFFInterpolationExportPython(slice2)
    chiExportFieldFunctionToVTK(fftemp,"ZPhi")
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZPFFI00.py")
end


