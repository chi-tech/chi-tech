if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/PolyFile.obj")

--############################################### Extract edges from surface mesh
loops,loop_count = chiSurfaceMeshGetEdgeLoops(newSurfMesh)

line_mesh = {};
line_mesh_count = 0;

for k=3,3 do
    split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
    for m=1,split_count do
        line_mesh_count = line_mesh_count + 1;
        line_mesh[line_mesh_count] = chiLineMeshCreateFromLoop(split_loops,m-1);
    end

end

--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
for k=1,line_mesh_count do
    chiRegionAddLineBoundary(region1,line_mesh[k]);
end

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_TRIANGLE);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);

chiSurfaceMesherSetProperty(MAX_AREA,1/20/20/4.0)
chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)


NZ=2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.5,NZ,"Charlie");
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.5,NZ,"Charlie");

--chiVolumeMesherSetProperty(PARTITION_Z,4)
--chiVolumeMesherSetProperty(PARTITION_Z,1)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--chiRegionExportMeshToPython(region1,"YMeshTT.py",true)

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");
materials[1] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)


chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[1],SCALAR_VALUE,SINGLE_VALUE,10000)


--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
fftemp = chiSolverAddFieldFunction(phys1,"Temperature")
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-4)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,2,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,3,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,4,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,5,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,6,DIFFUSION_DIRICHLET,0.0)


--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Set derived geometry
slice1 = chiFFInterpolationCreate(SLICE)
--chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.008,0.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0.0,0.0,1.0)
chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,0.0,-1.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,1.0,0.0,0.0)
chiFFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fftemp)

slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp)

line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-0.999,0.0,0.5)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT, 0.999,0.0,0.5)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(slice1)
chiFFInterpolationExecute(slice1)
chiFFInterpolationExportPython(slice1)

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)
chiFFInterpolationExportPython(slice2)

chiFFInterpolationInitialize(line0)
chiFFInterpolationExecute(line0)
chiFFInterpolationExportPython(line0)


if (chi_location_id == 0) then
    --local handle = io.popen("python ZPFFI10.py")
    local handle = io.popen("python ZLFFI20.py")
    print("Execution completed")
end