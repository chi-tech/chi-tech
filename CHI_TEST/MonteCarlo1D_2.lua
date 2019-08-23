chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)

--############################################### Setup mesh
tmesh = chiMeshHandlerCreate()

mesh={}
N=500
L=30.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
line_mesh = chiLineMeshCreateFromArray(mesh)


region0 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);


--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

chiVolumeMesherSetProperty(PARTITION_Z,6)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)




--############################################### Setup mesh
mmesh = chiMeshHandlerCreate()

mesh={}
N=50
L=30.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
line_mesh = chiLineMeshCreateFromArray(mesh)


region1 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);


--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

chiVolumeMesherSetProperty(MESH_GLOBAL,true)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)





--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_3_170.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,PDT_XSFILE,"xs_3_170.data")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)


src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics
chiMeshHandlerSetCurrent(mmesh)
phys0 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys0,region1)

chiMonteCarlonCreateSource(phys0,MC_BNDRY_SRC,1);

chiMonteCarlonSetProperty(phys0,MC_NUM_PARTICLES,1e6)
chiMonteCarlonSetProperty(phys0,MC_TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys0,MC_TALLY_MERGE_INTVL,100e3)
chiMonteCarlonSetProperty(phys0,MC_SCATTERING_ORDER,10)
chiMonteCarlonSetProperty(phys0,MC_MONOENERGETIC,false)
chiMonteCarlonSetProperty(phys0,MC_FORCE_ISOTROPIC,false)
chiMonteCarlonSetProperty(phys0,MC_TALLY_MULTIPLICATION_FACTOR,0.5)

chiMonteCarlonInitialize(phys0)
chiMonteCarlonExecute(phys0)


--############################################### Setup Physics
chiMeshHandlerSetCurrent(tmesh)
phys1 = chiNPTransportCreateSolver()
chiSolverAddRegion(phys1,region0)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiNPTCreateGroup(phys1)
end

--========== ProdQuad
pqaud = chiCreateProductQuadrature(GAUSS_LEGENDRE,40)

--========== Groupset def
gs0 = chiNPTCreateGroupset(phys1)
chiNPTGroupsetAddGroups(phys1,gs0,0,num_groups-1)
chiNPTGroupsetSetQuadrature(phys1,gs0,pquad)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
chiNPTSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiNPTSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiNPTSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiNPTSetProperty(phys1,SCATTERING_ORDER,5)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs0,NPT_CLASSICRICHARDSON)
--chiNPTSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs1,NPT_CLASSICRICHARDSON)
chiNPTSetProperty(phys1,GROUPSET_TOLERANCE,gs0,1.0e-6)
--chiNPTSetProperty(phys1,GROUPSET_MAXITERATIONS,gs0,3)
chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs0,100)
chiNPTSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs1,100)
chiNPTSetProperty(phys1,SWEEP_EAGER_LIMIT,62000)

chiNPTSetProperty(phys1,GROUPSET_SUBSETS,gs0,5)
chiNPTSetProperty(phys1,GROUPSET_SUBSETS,gs1,7)

chiNPTInitialize(phys1)
chiNPTExecute(phys1)






fflist,count = chiNPTGetScalarFieldFunctionList(phys1)

--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, 30.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

for k=1,2 do
    chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,k-1)
end

for k=1,2 do
    chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[k])
end


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)


--


if (chi_location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end
