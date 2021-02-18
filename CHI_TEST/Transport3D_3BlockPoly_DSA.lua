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

chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

NZ=1
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--10.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--20.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--30.0
chiVolumeMesherSetProperty(EXTRUSION_LAYER,10.0,NZ,"Charlie");--40.0

chiVolumeMesherSetProperty(PARTITION_Z,1);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);
chiVolumeMesherSetProperty(MESH_GLOBAL,false);
chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--chiRegionExportMeshToPython(region1,
--        "YMesh"..string.format("%d",chi_location_id)..".py",false)

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-10.0,10.0,-10.0,10.0,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

--chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_air50RH.data")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

--chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)




--############################################### Setup Physics
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,false)
pquad1 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8,false)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)

cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,62)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,5)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2)
-- chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")


gs1 = chiLBSCreateGroupset(phys1)

cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,63,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,5)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2,false," ")
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-6,false," ")

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
--bsrc[1] = 1.0
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMAX,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,INCIDENT_ISOTROPIC,bsrc);



chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi","Phi")
end
