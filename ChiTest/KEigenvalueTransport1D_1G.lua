-- 1D KEigen solver test with Vacuum BC.
-- SDM: PWLD
-- Test: Final k-eigenvalue: 0.997501
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

chiMPIBarrier()

-- ##################################################
-- ##### Parameters #####
-- ##################################################

-- Mesh variables
if (L == nil) then L = 100.0 end
if (n_cells == nil) then n_cells = 50 end

-- Transport angle information
if (n_angles == nil) then n_angles = 16 end 
if (scat_order == nil) then scat_order = 0 end

-- k-eigenvalue iteration parameters
if (max_k_iters == nil) then max_k_iters = 5000 end
if (k_tol == nil) then k_tol = 1e-8 end

-- Source iteration parameters
if (max_si_iters == nil) then max_si_iters = 500 end
if (si_tol == nil) then si_tol = 1e-4 end

-- Delayed neutrons
if (use_precursors == nil) then use_precursors = true end

-- Total cross section
if (sigma_t == nil) then sigma_t = 1.0 end

-- NOTE: For command line inputs, specify as:
--       variable=[[argument]]



-- ##################################################
-- ##### Run problem #####
-- ##################################################

--############################################### Setup mesh
chiMeshHandlerCreate()
nodes = {}
dx = L/n_cells
for i=0,n_cells do
  nodes[i+1] = i*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
chiVolumeMesherExecute()

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Fissile Material")

chiPhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                              CHI_XSFILE, "ChiTest/simple_fissile.cxs")

--############################################### Setup Physics
-- Define solver
phys = chiKEigenvalueLBSCreateSolver()

-- Add region and discretization
chiSolverAddRegion(phys, region)
chiLBSSetProperty(phys,DISCRETIZATION_METHOD,PWLD)

-- Create quadrature and define scattering order
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,n_angles)
chiLBSSetProperty(phys,SCATTERING_ORDER,scat_order)

-- Create groups
for g=0, num_groups - 1 do
    chiLBSCreateGroup(phys)
end

-- Create groupset
gs = chiLBSCreateGroupset(phys)
chiLBSGroupsetAddGroups(phys,gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys,gs,pquad)
chiLBSGroupsetSetMaxIterations(phys,gs,max_si_iters)
chiLBSGroupsetSetResidualTolerance(phys,gs,si_tol)
chiLBSGroupsetSetIterativeMethod(phys,gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetAngleAggregationType(phys,gs,LBSGroupset.ANGLE_AGG_SINGLE)

-- Additional parameters
chiLBSSetMaxKIterations(phys,max_k_iters)
chiLBSSetKTolerance(phys,k_tol)
chiLBSSetUsePrecursors(phys,use_precursors)
chiLBSSetProperty(phys,VERBOSE_INNER_ITERATIONS,false)
chiLBSSetProperty(phys,VERBOSE_OUTER_ITERATIONS,false)

--############################################### Initialize and Execute Solver
chiKEigenvalueLBSInitialize(phys)
chiKEigenvalueLBSExecute(phys)

--############################################### Get field functions
--############################################### Line plot
--############################################### Volume integrations
--############################################### Exports
--############################################### Plots