-- 2D Transport test with point source Multigroup Adjoint Response
-- SDM: PWLD
-- Test: Inner-product=3.30607e-06
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
mesh,region0 = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)
chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,0.0,0.8*L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)



----############################################### Set Material IDs
vol0b = chiLogicalVolumeCreate(RPP,-0.166666+2.5,0.166666+2.5,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0b,0)

vol1b = chiLogicalVolumeCreate(RPP,-1+2.5,1+2.5,0.9*L,L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1b,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 10
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.1*20,0.8)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0

--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 2)
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

--############################################### Add point source
chiLBSAddPointSource(phys1, 1.25 - 0.5*ds, 1.5*ds, 0.0, src)

ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

chiSolverInitialize(ss_solver)
--chiSolverExecute(ss_solver)

chiLBSReadFluxMoments(phys1, "Adjoint2D_3b_adjoint")
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
    os.execute("rm Adjoint2D_3b_adjoint*")
end
