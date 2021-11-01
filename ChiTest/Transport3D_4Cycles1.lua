-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=3.74343e-04
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

unpart_mesh = chiUnpartitionedMeshFromWavefrontOBJ(
        "ChiResources/TestObjects/Square2x2_partition_cyclic3.obj")

region1 = chiRegionCreate()

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER,
                      ExtruderTemplateType.UNPARTITIONED_MESH,
                      unpart_mesh);

NZ=2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.4
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.8
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.6

chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)
chiVolumeMesherSetKBAPartitioningPxPyPz(2,2,1)
chiVolumeMesherSetKBACutsX({0.0})
chiVolumeMesherSetKBACutsY({0.0})

chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 21
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_graphite_pure.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_graphite_pure.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics

phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,6, 6)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,20)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--
--gs1 = chiLBSCreateGroupset(phys1)
--chiLBSGroupsetAddGroups(phys1,gs1,63,167)
--chiLBSGroupsetSetQuadrature(phys1,gs1,pquad)

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,LBSBoundaryTypes.REFLECTING);

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
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[20])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")

    line = chiFFInterpolationCreate(LINE)
    chiFFInterpolationSetProperty(line,LINE_FIRSTPOINT,0.0,-1.0,0.5)
    chiFFInterpolationSetProperty(line,LINE_SECONDPOINT,0.0, 1.0,0.5)
    chiFFInterpolationSetProperty(line,LINE_NUMBEROFPOINTS,1000)
    chiFFInterpolationSetProperty(line,ADD_FIELDFUNCTION,fflist[2])

    chiFFInterpolationInitialize(line)
    chiFFInterpolationExecute(line)

    chiFFInterpolationExportPython(line,"Line")
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end

