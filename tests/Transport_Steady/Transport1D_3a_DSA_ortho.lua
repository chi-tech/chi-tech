-- 1D LinearBSolver test of a block of graphite with an air cavity. DSA and TG
-- SDM: PWLD
-- Test: WGS groups [0-62] Iteration    28 Residual 6.74299e-07 CONVERGED
-- and   WGS groups [63-167] Iteration    39 Residual 8.73816e-07 CONVERGED
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=1000
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

--chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

vol1 = chiLogicalVolumeCreate(RPP,-10.0,10.0,-10.0,10.0,-10.0,10.0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


--num_groups = 1
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        SIMPLEXS1,num_groups,1.0,0.999)
num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"tests/Transport_Steady/xs_graphite_pure.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"tests/Transport_Steady/xs_air50RH.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateSolver()

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
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
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON_CYCLES)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)

chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2,false,"")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,"")


gs1 = chiLBSCreateGroupset(phys1)

cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,63,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON_CYCLES)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)

chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-2,false,"")
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,"")

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
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Exports
if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi","Phi")
end

--############################################### Plots

--############################################### Stats N=50
-- Case 0a 1MPI GMRES(30)
-- GS0  53 its GS1 240 its

-- Case 1a 1MPI GMRES(30) WGDSA0 WGDSA1
-- GS0  27 its GS1  62 its

-- Case 2a 1MPI GMRES(30) WGDSA0 WGDSA1 TG1
-- GS0  27 its GS1  33 its


-- Case 0a 1MPI GMRES(60)
-- GS0  36 its GS1 70 its

-- Case 1a 1MPI GMRES(60) WGDSA0 WGDSA1
-- GS0  27 its GS1  42 its

-- Case 2a 1MPI GMRES(60) WGDSA0 WGDSA1 TG1
-- GS0  27 its GS1  32 its


-- GMRES(10)
-- Case 0b 1MPI GMRES(10)
-- GS0  96 its GS1 357 its

-- Case 1b 1MPI GMRES(10) WGDSA0 WGDSA1
-- GS0  116 its GS1  inf its

-- Case 2b 1MPI GMRES(10) WGDSA0 WGDSA1 TG1
-- GS0  116 its GS1  86 its

-- 8 MPI

-- GMRES(60)
-- Case 0c 8MPI GMRES(60)
-- GS0  39 its GS1 117 its

-- Case 1c 8MPI GMRES(60) WGDSA0 WGDSA1
-- GS0  32 its GS1 60 its

-- Case 2c 8MPI GMRES(60) WGDSA0 WGDSA1 TG1
-- GS0  32 its GS1  41 its


-- GMRES(30)
-- Case 0c 8MPI GMRES(30)
-- GS0  53 its GS1 239 its

-- Case 1c 8MPI GMRES(30) WGDSA0 WGDSA1
-- GS0  47 its GS1  110 its

-- Case 2c 8MPI GMRES(30) WGDSA0 WGDSA1 TG1
-- GS0  47 its GS1  56 its


-- GMRES(20)
-- Case 0c 8MPI GMRES(20)
-- GS0 109 its GS1 300 its

-- Case 1c 8MPI GMRES(20) WGDSA0 WGDSA1
-- GS0  57 its GS1 155 its

-- Case 2c 8MPI GMRES(20) WGDSA0 WGDSA1 TG1
-- GS0  57 its GS1  62 its

-- GMRES(10)
-- Case 0d 8MPI GMRES(10)
-- GS0  92 its GS1 410 its

-- Case 1d 8MPI GMRES(10) WGDSA0 WGDSA1
-- GS0  92 its GS1  inf its

-- Case 2d 8MPI GMRES(10) WGDSA0 WGDSA1 TG1
-- GS0  92 its GS1  98 its



--=============================================== BICGSTAB comparisons
-- Case 0a 1MPI GMRES(10)
-- GS0  87 its GS1 262 its TTS=56

-- Case 0a 1MPI GMRES(30)
-- GS0  54 its GS1 244 its TTS=46

-- Case 0a 1MPI GMRES(60)
-- GS0  37 its GS1 102 its TTS=21

-- BICGSTAB
-- Case 0a 1MPI BICGSTAB
-- GS0  34 its GS1 73 its TTS=30