print("############################################### LuaTest")
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

--############################################### Extract edges from surface mesh
-- loops,loop_count = chiSurfaceMeshGetEdgeLoopsPoly(newSurfMesh)
--
-- line_mesh = {};
-- line_mesh_count = 0;
--
-- for k=1,loop_count do
--     split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
--     for m=1,split_count do
--         line_mesh_count = line_mesh_count + 1;
--         line_mesh[line_mesh_count] =
--         chiLineMeshCreateFromLoop(split_loops,m-1);
--     end
--
-- end

--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
-- chiRegionAddEmptyBoundary(region1);
-- for k=1,line_mesh_count do
--     chiRegionAddLineBoundary(region1,line_mesh[k]);
-- end

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);

chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,-0.5)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.99999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

e_bndry = chiRegionAddEmptyBoundary(region1);
w_bndry = chiRegionAddEmptyBoundary(region1);
n_bndry = chiRegionAddEmptyBoundary(region1);
s_bndry = chiRegionAddEmptyBoundary(region1);

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)



--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
fftemp = chiSolverAddFieldFunction(phys1,"Temperature")
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLD_MIP);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-8)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,e_bndry,DIFFUSION_VACUUM)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,w_bndry,DIFFUSION_VACUUM)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,n_bndry,DIFFUSION_REFLECTING)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,s_bndry,DIFFUSION_REFLECTING)

--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)


ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

if (master_export == nil) then
    chiFFInterpolationExportPython(slice2)
    chiExportFieldFunctionToVTK(fftemp,"ZPhi")
end

if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZPFFI00.py")
end


