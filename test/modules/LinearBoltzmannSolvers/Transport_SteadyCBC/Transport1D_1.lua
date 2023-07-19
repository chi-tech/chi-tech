-- 1D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.49903 and 7.18243e-4
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
  chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=20
L=10.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
  k=i-1
  mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
RPP = chi_mesh.RPPLogicalVolume
lv0 = RPP.Create({infx=true, infy=true, zmin = -1.0, zmax=L/2})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, lv0, 1)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_3_170.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_3_170.cxs")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
  src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE,1)
lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 62},
      --groups_from_to = {0, num_groups-1},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
    {
      groups_from_to = {63, num_groups-1},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  options =
  {
    scattering_order = 5,
    save_angular_flux = true
  },
  sweep_type = "CBC"
  --sweep_type = "AAH"
}

--[0]  00:00:00 Computing b
--[0]  00:00:00 WGS groups [0-0] Iteration     0 Residual         1
--[0]  00:00:00 WGS groups [0-0] Iteration     1 Residual 0.0570112
--[0]  00:00:00 WGS groups [0-0] Iteration     2 Residual 0.00135648
--[0]  00:00:00 WGS groups [0-0] Iteration     3 Residual 2.18501e-05
--[0]  00:00:00 WGS groups [0-0] Iteration     4 Residual 9.43009e-07 CONVERGED


--[0]  Number of lagged angular unknowns: 0(0%)
--[0]  00:00:00 Computing b
--[0]  00:00:00 WGS groups [0-62] Iteration     0 Residual         1
--[0]  00:00:00 WGS groups [0-62] Iteration     1 Residual  0.117618
--[0]  00:00:00 WGS groups [0-62] Iteration     2 Residual 0.0182073
--[0]  00:00:00 WGS groups [0-62] Iteration     3 Residual 0.00599208
--[0]  00:00:00 WGS groups [0-62] Iteration     4 Residual 0.00238186
--[0]  00:00:00 WGS groups [0-62] Iteration     5 Residual 0.00100967
--[0]  00:00:00 WGS groups [0-62] Iteration     6 Residual 0.000451957
--[0]  00:00:00 WGS groups [0-62] Iteration     7 Residual 0.000208153
--[0]  00:00:00 WGS groups [0-62] Iteration     8 Residual 9.6763e-05
--[0]  00:00:00 WGS groups [0-62] Iteration     9 Residual 4.63033e-05
--[0]  00:00:00 WGS groups [0-62] Iteration    10 Residual 2.09523e-05
--[0]  00:00:00 WGS groups [0-62] Iteration    11 Residual 1.04169e-05
--[0]  00:00:00 WGS groups [0-62] Iteration    12 Residual 4.94007e-06
--[0]  00:00:00 WGS groups [0-62] Iteration    13 Residual 2.2935e-06
--[0]  00:00:00 WGS groups [0-62] Iteration    14 Residual 1.15192e-06
--[0]  00:00:00 WGS groups [0-62] Iteration    15 Residual 5.29544e-07 CONVERGED
--[0]
--[0]
--[0]
--[0]          Set Src Time/sweep (s):        0.000273706
--[0]          Average sweep time (s):        7.85294e-05
--[0]          Chunk-Overhead-Ratio  :        0.142322
--[0]          Sweep Time/Unknown (ns):       31.1625
--[0]          Number of unknowns per sweep:  5040
--[0]
--[0]
--[0]
--[0]
--[0]  ********** Solving groupset 1 with KRYLOV_GMRES.
--[0]
--[0]  Quadrature number of angles: 2
--[0]  Groups 63 167
--[0]
--[0]  Total number of angular unknowns: 8400
--[0]  Number of lagged angular unknowns: 0(0%)
--[0]  00:00:00 Computing b
--[0]  00:00:00 WGS groups [63-167] Iteration     0 Residual         1
--[0]  00:00:00 WGS groups [63-167] Iteration     1 Residual  0.414602
--[0]  00:00:00 WGS groups [63-167] Iteration     2 Residual  0.185221
--[0]  00:00:00 WGS groups [63-167] Iteration     3 Residual 0.0939791
--[0]  00:00:00 WGS groups [63-167] Iteration     4 Residual 0.0452038
--[0]  00:00:00 WGS groups [63-167] Iteration     5 Residual  0.021578
--[0]  00:00:00 WGS groups [63-167] Iteration     6 Residual 0.0105877
--[0]  00:00:00 WGS groups [63-167] Iteration     7 Residual 0.00517366
--[0]  00:00:00 WGS groups [63-167] Iteration     8 Residual 0.00219249
--[0]  00:00:00 WGS groups [63-167] Iteration     9 Residual 0.000921049
--[0]  00:00:00 WGS groups [63-167] Iteration    10 Residual  0.000485
--[0]  00:00:00 WGS groups [63-167] Iteration    11 Residual 0.000297586
--[0]  00:00:00 WGS groups [63-167] Iteration    12 Residual 0.000174054
--[0]  00:00:00 WGS groups [63-167] Iteration    13 Residual 0.000102658
--[0]  00:00:00 WGS groups [63-167] Iteration    14 Residual 6.11625e-05
--[0]  00:00:00 WGS groups [63-167] Iteration    15 Residual 3.46333e-05
--[0]  00:00:00 WGS groups [63-167] Iteration    16 Residual 2.0036e-05
--[0]  00:00:00 WGS groups [63-167] Iteration    17 Residual 1.08396e-05
--[0]  00:00:00 WGS groups [63-167] Iteration    18 Residual 5.58818e-06
--[0]  00:00:00 WGS groups [63-167] Iteration    19 Residual 2.57691e-06
--[0]  00:00:00 WGS groups [63-167] Iteration    20 Residual 1.23802e-06
--[0]  00:00:00 WGS groups [63-167] Iteration    21 Residual 6.00498e-07 CONVERGED
--[0]
--[0]
--[0]
--[0]          Set Src Time/sweep (s):        0.0005209
--[0]          Average sweep time (s):        0.000138478
--[0]          Chunk-Overhead-Ratio  :        0.219152
--[0]          Sweep Time/Unknown (ns):       32.971
--[0]          Number of unknowns per sweep:  8400
--[0]
--[0]
--[0]  Max-value1=1.99518
--[0]  Max-value2=1.30704e-14



phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)


--############################################### Volume integrations
vol0 = chi_mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[100])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if (master_export == nil) then
  --chiFFInterpolationExportPython(cline)
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
  --local handle = io.popen("python3 ZLFFI00.py")
end
