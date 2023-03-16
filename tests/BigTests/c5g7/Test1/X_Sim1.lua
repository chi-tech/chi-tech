dofile("mesh/gmesh_coarse.lua")

--chiMeshHandlerExportMeshToVTK("ZMesh")
--os.exit()

dofile("materials/materials.lua")

--############################################### Setup Physics
phys1 = chiLBKESCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
--pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,1, 1)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1, cur_gs, LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-8)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,50)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,50)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-8,false)
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-8,false)


--############################################### Set boundary conditions
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.REFLECTING);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 50)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)
chiLBKESSetProperty(phys1, "K_EIGEN_METHOD", "nonlinear")

chiLBSSetProperty(phys1, USE_PRECURSORS, false)

chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)

-- chiLBSSetProperty(phys1, SWEEP_EAGER_LIMIT, 1e9)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

chiExportMultiFieldFunctionToVTK(fflist,"solutions/ZPhi")
