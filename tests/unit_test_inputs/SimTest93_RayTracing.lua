--############################################### Setup mesh
chiMeshHandlerCreate()

if (nmesh==nil) then nmesh = 11 end

mesh={}
N=nmesh
L=11
xmin = -L/2
--xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

--chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

chiVolumeMesherSetupOrthogonalBoundaries()

chiSimTest93_RayTracing()


--###############################################
--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
        TRANSPORT_XSECTIONS,
        SIMPLEXS0,1,0.27)



--############################################### Setup Physics
solver_name = "LBS"
phys1 = chiLBSCreateSolver(solver_name)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12*2*4, 12*4)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,0)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions

--############################################### Add point source
src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiLBSAddPointSource(phys1, 0.0, 0.0, 0.0, src)

--############################################### Set solver properties
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

ff_m0 = chiGetFieldFunctionHandleByName(solver_name.."_Flux_g000_m00")

chiExportMultiFieldFunctionToVTK({ff_m0},"SimTest_93_LBS_"..solver_name)
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm SimTest_93*")
end