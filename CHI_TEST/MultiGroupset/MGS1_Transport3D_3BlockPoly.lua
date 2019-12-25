chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2QuadsBlock.obj",true)

--############################################### Extract edges from surface mesh
loops,loop_count = chiSurfaceMeshGetEdgeLoopsPoly(newSurfMesh)

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
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);

chiSurfaceMesherSetProperty(MAX_AREA,1/20/20)
chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

NZ=2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--10.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--20.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--30.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--40.0

chiVolumeMesherSetProperty(PARTITION_Z,1);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);
chiVolumeMesherSetProperty(MESH_GLOBAL,false);

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

chiRegionExportMeshToPython(region1,
        "YMesh"..string.format("%d",chi_location_id)..".py",false)

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE)
chiPhysicsMaterialAddProperty(materials[2],SCALAR_VALUE)

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

--chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
--chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


chiPhysicsMaterialSetProperty(materials[1],SCALAR_VALUE,SINGLE_VALUE,13.7)
chiPhysicsMaterialSetProperty(materials[2],SCALAR_VALUE,SINGLE_VALUE,0.5)

num_groups = 168
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.01)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.01)
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS1,num_groups,0.02,0.5)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS1,num_groups,0.04,0.5)
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_graphite_pure.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_graphite_pure.data")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

src[1] = 1.0
--chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
--chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics

phys1 = chiNPTransportCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiNPTCreateGroup(phys1)
end

--========== ProdQuad
pqaud = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,true)
pqaud2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,true)

--========== Groupset def
gs0 = chiNPTCreateGroupset(phys1)
chiNPTGroupsetAddGroups(phys1,gs0,0,62)
chiNPTGroupsetSetQuadrature(phys1,gs0,pqaud)

gs1 = chiNPTCreateGroupset(phys1)
chiNPTGroupsetAddGroups(phys1,gs1,63,num_groups-1)
chiNPTGroupsetSetQuadrature(phys1,gs1,pqaud2)



--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0;
end
bsrc[1] = 1.0/4.0/math.pi;
chiNPTSetProperty(phys1,BOUNDARY_CONDITION,XMIN,INCIDENT_ISOTROPIC,bsrc);
--chiNPTSetProperty(phys1,BOUNDARY_CONDITION,XMAX,INCIDENT_ISOTROPIC,bsrc);
--chiNPTSetProperty(phys1,BOUNDARY_CONDITION,YMIN,INCIDENT_ISOTROPIC,bsrc);
--chiNPTSetProperty(phys1,BOUNDARY_CONDITION,YMAX,INCIDENT_ISOTROPIC,bsrc);
--chiNPTSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);
--chiNPTSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiNPTSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiNPTSetProperty(phys1,DISCRETIZATION_METHOD,PWL3D)
chiNPTSetProperty(phys1,SCATTERING_ORDER,1)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs0,NPT_CLASSICRICHARDSON)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs1,NPT_CLASSICRICHARDSON)
chiNPTSetProperty(phys1,GROUPSET_TOLERANCE,gs0,1.0e-6)
chiNPTSetProperty(phys1,GROUPSET_MAXITERATIONS,gs0,200)
chiNPTSetProperty(phys1,GROUPSET_MAXITERATIONS,gs1,200)
chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs0,100)
chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs1,100)


chiNPTInitialize(phys1)
chiNPTExecute(phys1)



fflist,count = chiNPTGetScalarFieldFunctionList(phys1)
slices = {}
lines = {}
for k=1,count do
    slices[k] = chiFFInterpolationCreate(SLICE)
    chiFFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,20.01)
    chiFFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
    --chiFFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
    --chiFFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
    --chiFFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
    chiFFInterpolationInitialize(slices[k])
    chiFFInterpolationExecute(slices[k])
    chiFFInterpolationExportPython(slices[k])

    lines[k] = chiFFInterpolationCreate(LINE)
    chiFFInterpolationSetProperty(lines[k],LINE_FIRSTPOINT,-20.0,0.0,20.01)
    chiFFInterpolationSetProperty(lines[k],LINE_SECONDPOINT,20.0,0.0,20.01)
    chiFFInterpolationSetProperty(lines[k],LINE_NUMBEROFPOINTS, 100)
    chiFFInterpolationSetProperty(lines[k],ADD_FIELDFUNCTION,fflist[k])

    chiFFInterpolationInitialize(lines[k])
    chiFFInterpolationExecute(lines[k])
    chiFFInterpolationExportPython(lines[k])
end

volumes = {}
for k=1,count do
    volumes[k] = chiFFInterpolationCreate(VOLUME)
    chiFFInterpolationSetProperty(volumes[k],ADD_FIELDFUNCTION,fflist[k])
    chiFFInterpolationInitialize(volumes[k])
    chiFFInterpolationExecute(volumes[k])
    value = chiFFInterpolationGetValue(volumes[k])

    if (chi_location_id == 0) then
        print("Value for group "..tostring(k-1).." "..tostring(value))
    end
end



if (chi_location_id == 0) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    local handle = io.popen("python ZPFFI90.py")
    local handle = io.popen("python ZLFFI100.py")

    print("Execution completed")
end