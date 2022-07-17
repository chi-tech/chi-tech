-- 2D Transport test with localized material source
-- SDM: PWLD
-- Test: QOI-value=5.04171e-08
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
mesh = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)
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
solver_name = "LBS"
phys1 = chiLBSCreateSolver(solver_name)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 2)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,500)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions

--############################################### Add point source
chiLBSAddPointSource(phys1, 1.25 - 0.5*ds, 1.5*ds, 0.0, src)

--############################################### Set solver properties
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

--############################################### Create QOIs
tvol0 = chiLogicalVolumeCreate(RPP,2.3333,2.6666,4.16666,4.33333,-1000,1000)
tvol1 = chiLogicalVolumeCreate(RPP,0.5   ,0.8333,4.16666,4.33333,-1000,1000)





--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)





--############################################### Get field functions
ff_m0 = chiGetFieldFunctionHandleByName(solver_name.."-Flux_g0_m0")
ff_m1 = chiGetFieldFunctionHandleByName(solver_name.."-Flux_g0_m1")
ff_m2 = chiGetFieldFunctionHandleByName(solver_name.."-Flux_g0_m2")


--############################################### Slice plot

--############################################### Volume integrations
QOI_value_sum = 0.0
for g=0,num_groups-1 do
    ff = chiGetFieldFunctionHandleByName(solver_name.."-Flux_g"..tostring(g).."_m0")
    ffi1 = chiFFInterpolationCreate(VOLUME)
    curffi = ffi1
    chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM)
    chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,tvol1)
    chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,ff)

    chiFFInterpolationInitialize(curffi)
    chiFFInterpolationExecute(curffi)
    QOI_value = chiFFInterpolationGetValue(curffi)

    chiLog(LOG_0,string.format("QOI-value["..tostring(g).."]= %.5e", QOI_value))

    QOI_value_sum = QOI_value_sum + QOI_value
end
chiLog(LOG_0,string.format("QOI-value[sum]= %.5e", QOI_value_sum))

--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_"..solver_name)
    chiExportFieldFunctionToVTKG(ff_m0, "ZPhi_"..solver_name)
end

--[0]  QOI-vallue[0]= 1.95637e-09
--[0]  QOI-vallue[1]= 5.13773e-09
--[0]  QOI-vallue[2]= 6.82248e-09
--[0]  QOI-vallue[3]= 7.26518e-09
--[0]  QOI-vallue[4]= 6.76479e-09
--[0]  QOI-vallue[5]= 5.73749e-09
--[0]  QOI-vallue[6]= 2.68214e-09
--[0]  QOI-vallue[7]= 1.17145e-09
--[0]  QOI-vallue[8]= 5.31336e-10
--[0]  QOI-vallue[9]= 3.59573e-10
