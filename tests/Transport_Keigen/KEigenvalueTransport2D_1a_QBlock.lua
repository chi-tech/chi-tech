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
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,50)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,50)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-8,false)
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-8,false)


--############################################### Set boundary conditions
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.REFLECTING);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,2)

chiSolverSetBasicOption(phys1, "K_EIGEN_METHOD", "power")
chiSolverSetBasicOption(phys1, "PI_MAX_ITS", 1000)
chiSolverSetBasicOption(phys1, "PI_K_TOL", 1.0e-10)

chiLBSSetProperty(phys1, USE_PRECURSORS, false)

chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--chiExportMultiFieldFunctionToVTK(fflist,"tests/BigTests/QBlock/solutions/Flux")

-- Reference value k_eff = 0.5969127