-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.50758 and 2.52527e-04
num_procs = 4





--############################################### Check num_procs
-- if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
--     chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
--                       "Expected "..tostring(num_procs)..
--                       ". Pass check_num_procs=false to override if possible.")
--     os.exit(false)
-- end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=20
L=5
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
vol1 = chiLogicalVolumeCreate(RPP,-1000,0,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"tests/simple_scatter.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"tests/simple_scatter.cxs")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
fac=1
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4*fac, 3*fac)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)
--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,0)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-8)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

-- gs0 = chiLBSCreateGroupset(phys1)
-- cur_gs = gs0
-- chiLBSGroupsetAddGroups(phys1,cur_gs,1,62)
-- chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
-- chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
-- chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
-- chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
-- chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
-- chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
-- chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

-- gs1 = chiLBSCreateGroupset(phys1)
-- cur_gs = gs1
-- chiLBSGroupsetAddGroups(phys1,cur_gs,63,167)
-- chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
-- chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
-- chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
-- chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
-- chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
-- chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
-- chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions
-- bsrc={}
-- for g=1,num_groups do
--     bsrc[g] = 0.0
-- end
-- bsrc[1] = 1.0/4.0/math.pi
-- chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
--                         LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

chiLBSComputeBalance(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Slice plot
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if master_export == nil then
    chiFFInterpolationExportPython(slice2)
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZPFFI00.py")
end