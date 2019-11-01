chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)

--############################################### Setup Transport mesh
tmesh = chiMeshHandlerCreate()

mesh={}
N=20
L=5.0
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

--chiVolumeMesherSetProperty(PARTITION_Z,6)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)




--############################################### Setup Monte-carlo mesh
mmesh = chiMeshHandlerCreate()

mesh={}
--N=20
L=5.0
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


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS1,1,1.0,0.5)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS1,1,1.0,0.5)

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)


src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)






--############################################### Setup Transport Physics
chiMeshHandlerSetCurrent(tmesh)
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region0)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,1)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs0,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
--chiLBSSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs0,LBS_CLASSICRICHARDSON)
--chiLBSSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs1,NPT_CLASSICRICHARDSON)
chiLBSSetProperty(phys1,GROUPSET_TOLERANCE,gs0,1.0e-6)
--chiLBSSetProperty(phys1,GROUPSET_MAXITERATIONS,gs0,3)
chiLBSSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs0,100)
chiLBSSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs1,100)
chiLBSSetProperty(phys1,SWEEP_EAGER_LIMIT,62000)

chiLBSSetProperty(phys1,GROUPSET_SUBSETS,gs0,5)
chiLBSSetProperty(phys1,GROUPSET_SUBSETS,gs1,7)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist1,count = chiLBSGetScalarFieldFunctionList(phys1)


--############################################### Setup Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys0 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys0,region0)

--chiMonteCarlonCreateSource(phys0,MC_BNDRY_SRC,1);
chiMonteCarlonCreateSource(phys0,MC_RESID_SRC,fflist1[1]);

chiMonteCarlonSetProperty(phys0,MC_NUM_PARTICLES,10e6)
chiMonteCarlonSetProperty(phys0,MC_TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys0,MC_TALLY_MERGE_INTVL,100e3)
chiMonteCarlonSetProperty(phys0,MC_SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys0,MC_MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys0,MC_FORCE_ISOTROPIC,true)
chiMonteCarlonSetProperty(phys0,MC_TALLY_MULTIPLICATION_FACTOR,0.5*N)

chiMonteCarlonInitialize(phys0)
chiMonteCarlonExecute(phys0)

--############################################### Setup ref Monte Carlo Physics
chiMeshHandlerSetCurrent(mmesh)
phys2 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys2,region1)

chiMonteCarlonCreateSource(phys2,MC_BNDRY_SRC,1);
--chiMonteCarlonCreateSource(phys2,MC_RESID_SRC,fflist1[1]);

chiMonteCarlonSetProperty(phys2,MC_NUM_PARTICLES,1e6)
chiMonteCarlonSetProperty(phys2,MC_TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys2,MC_TALLY_MERGE_INTVL,100e3)
chiMonteCarlonSetProperty(phys2,MC_SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys2,MC_MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys2,MC_FORCE_ISOTROPIC,true)
chiMonteCarlonSetProperty(phys2,MC_TALLY_MULTIPLICATION_FACTOR,0.5)

chiMonteCarlonInitialize(phys2)
chiMonteCarlonExecute(phys2)


fflist0,count = chiGetFieldFunctionList(phys0)
fflist1,count = chiLBSGetScalarFieldFunctionList(phys1)
fflist2,count = chiGetFieldFunctionList(phys2) --Fine mesh MC

--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, 5.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist0[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist1[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1])

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)


--


if (chi_location_id == 0) then
    --local handle = io.popen("python ZLFFI00.py")
end
