-- 3D Transport test with Vacuum BCs.
-- SDM: PWLD
-- Test: Max-value=1.02943
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromEnsightGold("ChiResources/TestObjects/Sphere.case")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiVolumeMesherSetProperty(VOLUMEPARTITION_X,2)
chiVolumeMesherSetProperty(VOLUMEPARTITION_Y,2)
chiVolumeMesherSetProperty(CUTS_X,0.0)
chiVolumeMesherSetProperty(CUTS_Y,0.0)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

num_groups = 5
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_graphite_pure.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_graphite_pure.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1]=1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics

phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad  = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 8)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
if (master_export == nil) then
    --chiLBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
end
--chiLBSGroupsetSetMaxIterations(phys1,cur_gs,10)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)

--############################################### Initialize and Execute Solver
chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Slice plot
--slices = {}
--for k=1,count do
--    slices[k] = chiFFInterpolationCreate(SLICE)
--    chiFFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,0.8001)
--    chiFFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
--    --chiFFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
--    --chiFFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
--    --chiFFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
--    chiFFInterpolationInitialize(slices[k])
--    chiFFInterpolationExecute(slices[k])
--    chiFFInterpolationExportPython(slices[k])
--end

--############################################### Volume integrations
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5e", maxval))

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[2])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end


