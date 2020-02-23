if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);

NZ=10
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");

chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)



--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
--fftemp = chiSolverAddFieldFunction(phys1,"Temperature")
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLD_MIP);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-6)

--############################################### Set boundary conditions
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,5,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,6,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,0,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,2,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,3,DIFFUSION_VACUUM)


--############################################### Initialize Solver
chiDiffusionInitialize(phys1)



chiDiffusionExecute(phys1)
--
----############################################### Set derived geometry
fflist,count = chiGetFieldFunctionList(phys1)
--slice1 = chiFFInterpolationCreate(SLICE)
--
--chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.008,0.0,0.0)
--chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0.0,0.0,1.0)
--chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,0.0,-1.0,0.0)
--chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,1.0,0.0,0.0)
--chiFFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])
--
--slice2 = chiFFInterpolationCreate(SLICE)
--chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.024)
--chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])
--
chiMPIBarrier();
line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.0,0.025)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.0,0.025)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fflist[1])
--
--chiFFInterpolationInitialize(slice1)
--chiFFInterpolationExecute(slice1)
--
--chiFFInterpolationInitialize(slice2)
--chiFFInterpolationExecute(slice2)


print("Yip2")
chiMPIBarrier();
chiFFInterpolationInitialize(line0)
print("Yip2a")
chiMPIBarrier();
chiFFInterpolationExecute(line0)

print("Yip2b")
chiMPIBarrier();
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

print("Yip3")
chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

print("Yip3")
chiMPIBarrier();
if (master_export == nil) then
    chiFFInterpolationExportPython(slice1)
    chiFFInterpolationExportPython(slice2)
    chiFFInterpolationExportPython(line0)
end

chiMPIBarrier()

print("Yip4")
chiMPIBarrier();
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZPFFI10.py")
    local handle = io.popen("python ZLFFI20.py")
    print("Execution completed")
end

if (master_export == nil) then
    chiExportFieldFunctionToVTK(fflist[1],"ZPhi3D","Temperature")
end
