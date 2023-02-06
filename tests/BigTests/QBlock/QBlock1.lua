dofile("mesh.lua")
dofile("materials.lua") --num_groups assigned here

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
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-8)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,100)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-8,false)
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-8,false)


--############################################### Set boundary conditions
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.REFLECTING);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,2)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 50)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-10)

chiLBSSetProperty(phys1, USE_PRECURSORS, false)

chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

chiExportMultiFieldFunctionToVTK(fflist,"solutions/Flux")
