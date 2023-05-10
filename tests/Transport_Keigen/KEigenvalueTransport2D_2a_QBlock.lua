-- 2D 2G KEigenvalue::Solver test using Power Iteration
-- Test: Final k-eigenvalue: 0.5969127

dofile("tests/Transport_Keigen/QBlock_mesh.lua")
dofile("tests/Transport_Keigen/QBlock_materials.lua") --num_groups assigned here

--############################################### Setup Physics
--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 4)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

phys1 = lbs.DiscOrdSteadyStateSolver.Create
({
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, num_groups-1},
      angular_quadrature_handle = pquad,
      inner_linear_method = "gmres",
      l_max_its = 10,
      l_abs_tol = 1.0e-6
    }
  }
})

chiLBSSetOptions(phys1,
  {
    boundary_conditions = { { name = "xmin", type = "reflecting"},
                            { name = "ymin", type = "reflecting"} },
    spatial_discretization = "pwld",
    scattering_order = 2,
    save_angular_flux = true,

    use_precursors = false,

    verbose_inner_iterations = false,
    verbose_outer_iterations = true,

    power_field_function_on = true,
    power_default_kappa = 1.0,
    power_normalization = -1.0, --Disabled
  })

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
--k_solver = lbs.XXPowerIterationKEigenSCDSA.Create
--({
--  lbs_solver_handle = phys1,
--  diff_accel_sdm = "pwld",
--  accel_pi_verbose = false
--})
--chiSolverExecute(k_solver)

k_solver = lbs.XXNonLinearKEigen.Create
({
  lbs_solver_handle = phys1
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
  for k=5, 5+num_groups-1 do
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
    input_dimension = 4+count,
    output_dimension = 1,
    lua_function_name = "ComputePower"
  })
})

chiFieldOperationExecute(op)


fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
table.insert(fflist, ff_power)
table.insert(fflist, chiGetFieldFunctionHandleByName("power_generation"))

chiExportMultiFieldFunctionToVTK(fflist,"ZPhi2a")

-- Reference value k_eff = 0.5969127