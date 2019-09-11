chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/PolyFile2.obj")

--############################################### Extract edges from surface mesh
loops,loop_count = chiSurfaceMeshGetEdgeLoops(newSurfMesh)

line_mesh = {};
line_mesh_count = 0;

for k=1,loop_count do
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
--chiSurfaceMesherSetProperty(PARTITION_X,2)
--chiSurfaceMesherSetProperty(PARTITION_Y,2)
--chiSurfaceMesherSetProperty(CUT_X,0.0)
--chiSurfaceMesherSetProperty(CUT_Y,0.0)
chiSurfaceMesherSetProperty(PARTITION_X,3)
chiSurfaceMesherSetProperty(PARTITION_Y,1)
chiSurfaceMesherSetProperty(CUT_X,-0.33)
--chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_X,0.33)

NZ=2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.4
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.8
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.6

chiVolumeMesherSetProperty(PARTITION_Z,2);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);
chiVolumeMesherSetProperty(MESH_GLOBAL,false);

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

chiRegionExportMeshToPython(region1,
        "YMesh"..string.format("%d",chi_location_id)..".py",true)

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");
materials[1] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[0],THERMAL_CONDUCTIVITY)
chiPhysicsMaterialSetProperty(materials[0],THERMAL_CONDUCTIVITY,SINGLE_VALUE,13.7)

chiPhysicsMaterialAddProperty(materials[1],THERMAL_CONDUCTIVITY)
chiPhysicsMaterialSetProperty(materials[1],THERMAL_CONDUCTIVITY,SINGLE_VALUE,0.5)


--chiPhysicsMaterialSetProperty(materials[0],THERMAL_CONDUCTIVITY,FROM_TABLE,12.7)

chiPhysicsMaterialAddProperty(materials[0],TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[0],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")



--############################################### Setup Physics

phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
grp[0] = chiLBSCreateGroup(phys1)
grp[1] = chiLBSCreateGroup(phys1)

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs0,0,1)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)


chiLBSSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)

bsrc={}
for g=1,2 do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
--bsrc[1] = 1.0
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,INCIDENT_ISOTROPIC,bsrc);


chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
slices = {}
for k=1,count do
    slices[k] = chiFFInterpolationCreate(SLICE)
    chiFFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,0.025)
    chiFFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
    --chiFFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
    --chiFFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
    --chiFFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
    chiFFInterpolationInitialize(slices[k])
    chiFFInterpolationExecute(slices[k])
    chiFFInterpolationExportPython(slices[k])
end

if (chi_location_id == 0) then
    print("Execution completed")
    os.execute("python ZPFFI00.py")
    os.execute("python ZPFFI11.py")
end

if (chi_location_id == 0) then
    print("Execution completed")
end