-- 2D transport test in axialsymmetric cylindrical geometry with
-- vacuum boundary condition - multigroup with DSA.
-- SDM: PWLD
-- Test: Max-valueG1=1.00000, Max-valueG2=0.25000
num_procs = 4
--Structured mesh




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--------------------------------------------------------------------------------
--  mesh
--------------------------------------------------------------------------------
chiMeshHandlerCreate()
dim = 2
length = {1, 2, }
ncells = {50, 100, }
nodes = {}
for d = 1, dim do
  delta = length[d]/ncells[d]
  nodes[d] = {}
  for i = 0, ncells[d] do
    nodes[d][i+1] = i*delta
  end
end
surf_mesh = chiMeshCreateUnpartitioned2DOrthoMesh(nodes[1], nodes[2])
chiVolumeMesherSetProperty(PARTITION_TYPE, PARMETIS)
chiVolumeMesherExecute()

vol0 = chiLogicalVolumeCreate(RPP, 0, length[1], 0, length[2], 0, 0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--------------------------------------------------------------------------------
--  materials
--------------------------------------------------------------------------------
ngrp = 2
sigmat = 20.0
ratioc = 0.4
source = {}
source[1] = sigmat * (1 - 0.5*ratioc)
for g = 2, ngrp do
  source[g] = 0
end

material0 = chiPhysicsAddMaterial("Material_0");
chiPhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(material0,ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialSetProperty(material0, TRANSPORT_XSECTIONS,
                              SIMPLEXS1, ngrp, sigmat, ratioc)
chiPhysicsMaterialSetProperty(material0, ISOTROPIC_MG_SOURCE,
                              FROM_ARRAY, source)

--------------------------------------------------------------------------------
--  physics
--------------------------------------------------------------------------------
phys0 = chiLBSCurvilinearCreateSolver(LBSCurvilinear.CYLINDRICAL)

--  angular quadrature
pquad = chiCreateCylindricalProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 8)

--  groups
groups = {}
for g = 1, ngrp do
  groups[g] = chiLBSCreateGroup(phys0)
end

--  groupsets
gs0 = chiLBSCreateGroupset(phys0)
chiLBSGroupsetAddGroups(phys0, gs0, 0, ngrp-1)
chiLBSGroupsetSetQuadrature(phys0, gs0, pquad)
chiLBSGroupsetSetAngleAggregationType(phys0, gs0, LBSGroupset.ANGLE_AGG_AZIMUTHAL)
chiLBSGroupsetSetIterativeMethod(phys0, gs0, NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys0, gs0, 1.0e-12)
chiLBSGroupsetSetMaxIterations(phys0, gs0, 100)
chiLBSGroupsetSetGMRESRestartIntvl(phys0, gs0, 30)
petsc_options = ""
chiLBSGroupsetSetWGDSA(phys0, gs0, 50, 1.0e-09, false, petsc_options)
--chiLBSGroupsetSetTGDSA(phys0, gs0, 50, 1.0e-09, false, petsc_options)

--  spatial discretisation
chiLBSSetProperty(phys0, DISCRETIZATION_METHOD, PWLD)

--  scattering order
chiLBSSetProperty(phys0, SCATTERING_ORDER, 0)

--------------------------------------------------------------------------------
--  boundary conditions
--------------------------------------------------------------------------------
dirichlet_value = {}
for g = 1, ngrp do
  dirichlet_value[g] = 0
end
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  XMIN, LBSBoundaryTypes.REFLECTING)
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  XMAX, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  YMIN, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  YMAX, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)

--------------------------------------------------------------------------------
--  solvers
--------------------------------------------------------------------------------
chiSolverInitialize(phys0)
chiSolverExecute(phys0)

--------------------------------------------------------------------------------
--  output
--------------------------------------------------------------------------------
--  field functions
fflist, count = chiLBSGetScalarFieldFunctionList(phys0)
if master_export == nil then
  chiExportFieldFunctionToVTKG(fflist,
                               "Transport2DCyl2-scalar_flux",
                               "scalar_flux")
end

--  volume integrations - energy group 1
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi, OPERATION, OP_MAX)
chiFFInterpolationSetProperty(curffi, LOGICAL_VOLUME, vol0)
chiFFInterpolationSetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-valueG1=%.5f", maxval))

--  volume integrations - energy group 2
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[2])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-valueG2=%.5f", maxval))
