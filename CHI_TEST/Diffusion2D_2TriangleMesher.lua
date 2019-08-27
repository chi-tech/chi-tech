print("############################################### LuaTest")
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

if (chi_location_id == 0) then
    tempSurfMesh = chiSurfaceMeshCreate();
    chiSurfaceMeshImportFromOBJFile(tempSurfMesh,
            "CHI_RESOURCES/TestObjects/PolyFile.obj")
    chiSurfaceMeshExportPolyFile(tempSurfMesh,"ZTestPoly.poly")
    command = "CHI_RESOURCES/Dependencies/triangle/triangle -pqa0.000625 "
    command = command .. "ZTestPoly.poly"
    os.execute(command)

    --chiSurfaceMesherExportToObj(meshedSurfMesh,"ZRemeshedSurface.obj")
end
chiMPIBarrier()
meshedSurfMesh = chiSurfaceMeshCreate();
print("Creating surface mersh")
chiSurfaceMeshImportFromTriangleFiles(meshedSurfMesh,"ZTestPoly")
--newSurfMesh = chiSurfaceMeshCreate();
--chiSurfaceMeshImportFromOBJFile(newSurfMesh,
--        "CHI_RESOURCES/TestObjects/PolyFile.obj")
newSurfMesh = meshedSurfMesh


--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);

chiSurfaceMesherSetProperty(MAX_AREA,1/20/20/4)
chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],THERMAL_CONDUCTIVITY)
chiPhysicsMaterialSetProperty(materials[0],THERMAL_CONDUCTIVITY,SINGLE_VALUE,1.0)



--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
fftemp = chiSolverAddFieldFunction(phys1,"Temperature")
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-4)

--############################################### Set boundary conditions
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIRICHLET,0.1)

--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)
chiFFInterpolationExportPython(slice2)

line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.0,0.024)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.0,0.024)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(line0)
chiFFInterpolationExecute(line0)
chiFFInterpolationExportPython(line0)

if (chi_location_id == 0) then
    local handle = io.popen("python ZPFFI00.py")
end

