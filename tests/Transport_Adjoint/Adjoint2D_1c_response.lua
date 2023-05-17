-- 2D Transport test with localized material source Adjoint response
-- SDM: PWLD
-- Test: Inner-product=1.38405e-05
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
tmesh = chiMeshHandlerCreate()

nodes={}
N=60
L=5.0
ds=L/N
xmin=0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end
mesh = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)
chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,0.0,0.8*L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)



----############################################### Set Material IDs
vol0b = chiLogicalVolumeCreate(RPP,-0.166666+2.5,0.166666+2.5,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0b,0)

vol2 = chiLogicalVolumeCreate(RPP,-0.166666+2.5,0.166666+2.5,0.0,2*0.166666,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol2,2)

vol1b = chiLogicalVolumeCreate(RPP,-1+2.5,1+2.5,0.9*L,L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1b,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");
materials[3] = chiPhysicsAddMaterial("Test Material3");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[3],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[3],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.1*20,0.8)
chiPhysicsMaterialSetProperty(materials[3],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.3*20,0.0)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 3.0
chiPhysicsMaterialSetProperty(materials[3],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)


--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,48, 6)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud0, 4.0*math.pi)

lbs_block =
{
    num_groups = num_groups,
    groupsets =
    {
        {
            groups_from_to = {0, num_groups-1},
            angular_quadrature_handle = pquad0,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    }
}

lbs_options =
{
    scattering_order = 1,
}

--############################################### Initialize and Execute Solver
phys1 = lbs.DiscreteOrdinatesAdjointSolver.Create(lbs_block)
chiLBSSetOptions(phys1, lbs_options)

--############################################### Create QOIs
tvol0 = chiLogicalVolumeCreate(RPP,2.3333,2.6666,4.16666,4.33333,-1000,1000)
tvol1 = chiLogicalVolumeCreate(RPP,0.5   ,0.8333,4.16666,4.33333,-1000,1000)

chiAdjointSolverAddResponseFunction(phys1,"QOI0",tvol0)
chiAdjointSolverAddResponseFunction(phys1,"QOI1",tvol1)
chiSolverSetBasicOption(phys1, "REFERENCE_RF", "QOI1")

ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

chiSolverInitialize(ss_solver)
--chiSolverExecute(ss_solver)

chiLBSReadFluxMoments(phys1, "Adjoint2D_1b_adjoint")
value = chiAdjointSolverComputeInnerProduct(phys1)
chiLog(LOG_0,string.format("Inner-product=%.5e", value))

--############################################### Get field functions
ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")


--############################################### Slice plot

--############################################### Volume integrations



--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_LBAdjointResponse")
end


--############################################### Cleanup
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm Adjoint2D_1b_adjoint*")
end
