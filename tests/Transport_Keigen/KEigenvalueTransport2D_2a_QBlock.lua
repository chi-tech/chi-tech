-- 2D 2G KEigenvalue::Solver test using Power Iteration
-- Test: Final k-eigenvalue: 0.5969127

dofile("tests/Transport_Keigen/QBlock_mesh.lua")
dofile("tests/Transport_Keigen/QBlock_materials.lua") --num_groups assigned here

--############################################### Setup Physics
phys1 = chiLBKESCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 4)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON_CYCLES)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-10)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,10)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,50)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-8,false)
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-8,false)


--############################################### Set boundary conditions
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.REFLECTING);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,2)
chiLBSSetProperty(phys1,SAVE_ANGULAR_FLUX,true)

chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "power")
chiSolverSetBasicOption(phys1, "PI_MAX_ITS", 1000)
chiSolverSetBasicOption(phys1, "PI_K_TOL", 1.0e-10)

chiLBSSetProperty(phys1, USE_PRECURSORS, false)

chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
k_solver = lbs.XXPowerIterationKEigenSCDSA.Create
({
    lbs_solver_handle = phys1,
    diff_accel_sdm = "pwlc",
    accel_pi_verbose = false
})
chiSolverExecute(k_solver)
--chiSolverExecute(phys1)

--SCDSA 40 * 4 0.597 7278
--SMM   40 * 4 0.606 3770

ff_power = chi_physics.FieldFunctionGridBased.Create
({
    name = "fission_power",
    sdm_type = "PWLD",
    initial_value = 0.0
})

function ComputePower(params)
    local x = params[1]
    local y = params[2]
    local z = params[3]

    local mat_id = tostring(math.floor(params[4]))

    local xs_vals = chiPhysicsTransportXSGet(xs[mat_id])

    if (not xs_vals.is_fissionable) then return {0.0} end

    local FR = 0.0
    local g=0
    for k=5, rawlen(params) do
        g = g + 1
        FR = FR + xs_vals.sigma_f[g] * params[k]
    end

    return {FR}
end

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

op = chi_physics.field_operations.MultiFieldOperation.Create
({
    result_field_handle = ff_power,
    dependent_field_handles = fflist,
    function_handle = chi_math.functions.LuaDimAToDimB.Create
    ({
        input_dimension = 4+num_groups,
        output_dimension = 1,
        lua_function_name = "ComputePower"
    })
})

chiFieldOperationExecute(op)


fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
table.insert(fflist, ff_power)

chiExportMultiFieldFunctionToVTK(fflist,"ZPhi2a")

-- Reference value k_eff = 0.5969127