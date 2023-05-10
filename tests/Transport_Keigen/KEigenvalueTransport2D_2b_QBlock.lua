-- 2D 2G KEigenvalue::Solver test using Power Iteration
-- Test: Final k-eigenvalue: 0.5969127

dofile("tests/Transport_Keigen/QBlock_mesh.lua")
dofile("tests/Transport_Keigen/QBlock_materials.lua") --num_groups assigned here

--############################################### Setup Physics
--========== ProdQuad
pquadB = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 4)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaudB, 4.0*math.pi)
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 4)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, num_groups-1},
      angular_quadrature_handle = pquadB,
      inner_linear_method = "gmres",
      l_max_its = 10,
      l_abs_tol = 1.0e-6
    }
  }
}

lbs_options =
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

  field_function_prefix_option = "solver_name"
}

lbs_block.name = "MIP"
phys0 = lbs.DiffusionDFEMSolver.Create(lbs_block)
lbs_block.name = "DO"
lbs_block.groupsets[1].angular_quadrature_handle = pquad
phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

chiLBSSetOptions(phys0, lbs_options)
chiLBSSetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys0)
chiSolverInitialize(phys1)

k_solver0 = lbs.XXPowerIterationKEigen.Create({ lbs_solver_handle = phys0, })
k_solver1 = lbs.XXPowerIterationKEigenSCDSA.Create({ lbs_solver_handle = phys1,
                                             reinit_phi_1=true })
--k_solver0 = lbs.XXNonLinearKEigen.Create({ lbs_solver_handle = phys0, })
--k_solver1 = lbs.XXNonLinearKEigen.Create({ lbs_solver_handle = phys1,
--                                           reinit_phi_1=true})
chiSolverExecute(k_solver0)

op = chi_physics.field_operations.FieldCopyOperation.Create
({
  from = chiGetFieldFunctionHandleByName("MIP_phi_g000_m00"),
  to = chiGetFieldFunctionHandleByName("DO_phi_g000_m00")
})
chiFieldOperationExecute(op)
chiLBSSetPhiFromFieldFunction(phys1,{m_ids={0}})

chiSolverExecute(k_solver1)

--k_solver = lbs.XXPowerIterationKEigenSCDSA.Create
--({
--  lbs_solver_handle = phys1,
--  diff_accel_sdm = "pwld",
--  accel_pi_verbose = false
--bnhp
--k_solver = lbs.XXNonLinearKEigen.Create
--({
--  lbs_solver_handle = phys1
--})
--chiSolverExecute(k_solver)

fflist0,count = chiLBSGetScalarFieldFunctionList(phys0)
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

table.insert(fflist0, chiGetFieldFunctionHandleByName("MIP_power_generation"))
table.insert(fflist, chiGetFieldFunctionHandleByName("DO_power_generation"))

chiExportMultiFieldFunctionToVTK(fflist,"ZPhi2b")
chiExportMultiFieldFunctionToVTK(fflist0,"ZPhi2bb")


-- Reference value k_eff = 0.5969127