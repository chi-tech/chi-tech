chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

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
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.4
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--0.8
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.2
chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2*NZ,NZ,"Charlie");--1.6

chiVolumeMesherSetProperty(PARTITION_Z,1);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);
chiVolumeMesherSetProperty(MESH_GLOBAL,false);

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--chiRegionExportMeshToPython(region1,
--        "YMesh"..string.format("%d",chi_location_id)..".py",false)

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
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_graphite_pure.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_graphite_pure.data")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics

phys1 = chiNPTransportCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiNPTCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,16, 16)

--========== Groupset def
gs0 = chiNPTCreateGroupset(phys1)
chiNPTGroupsetAddGroups(phys1,gs0,0,20)
chiNPTGroupsetSetQuadrature(phys1,gs0,pquad)
--
--gs1 = chiNPTCreateGroupset(phys1)
--chiNPTGroupsetAddGroups(phys1,gs1,63,167)
--chiNPTGroupsetSetQuadrature(phys1,gs1,pquad)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
chiNPTSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiNPTSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiNPTSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs0,NPT_CLASSICRICHARDSON)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs1,NPT_CLASSICRICHARDSON)
chiNPTSetProperty(phys1,GROUPSET_TOLERANCE,gs0,1.0e-6)
--chiNPTSetProperty(phys1,GROUPSET_MAXITERATIONS,gs0,3)
chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs0,100)
--chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs1,100)

chiNPTSetProperty(phys1,GROUPSET_SUBSETS,gs0,3)
--chiNPTSetProperty(phys1,GROUPSET_SUBSETS,gs1,5)

chiNPTInitialize(phys1)
chiNPTExecute(phys1)



fflist,count = chiNPTGetScalarFieldFunctionList(phys1)
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

if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end

if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")
end

