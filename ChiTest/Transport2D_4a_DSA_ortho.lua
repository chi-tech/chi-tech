-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Nonce
num_procs = 4





----############################################### Check num_procs
--if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
--    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
--                      "Expected "..tostring(num_procs)..
--                      ". Pass check_num_procs=false to override if possible.")
--    os.exit(false)
--end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=50
L=100
--N=10
--L=200e6
xmin = -L/2
--xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
--vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
--chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--vol1 = chiLogicalVolumeCreate(RPP,-10.0,10.0,-10.0,10.0,-1000,1000)
--chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
--materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
--chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
--chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


--num_groups = 1
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        SIMPLEXS1,num_groups,1.0,0.999)
num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_graphite_pure.cxs")
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
--        CHI_XSFILE,"ChiTest/xs_air50RH.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
--src[1] = 0.0
--chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateSolver()

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,false)
--pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,false)
pquad1 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8,false)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)

cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,62)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
if (nodsa == nil) then
    chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2)
end
-- chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")


gs1 = chiLBSCreateGroupset(phys1)

cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,63,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
if (nodsa == nil) then
    chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2,false," ")
    chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-6,false," ")
end

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
--bsrc[1] = 1.0
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMAX,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,INCIDENT_ISOTROPIC,bsrc);

--############################################### Initialize and Execute Solver
chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Exports
if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi","Phi")
end

--############################################### Plots

--############################################### Stats
-- 1MPI RICH
-- GS0 140 its GS1 inf its

-- 1MPI RICH WGDSA0 WGDSA1
-- GS0  42 its GS1 463 its

-- 1MPI RICH WGDSA0        TG1
-- GS0  42 its GS1 --- its