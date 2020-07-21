chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



---############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
meshz={}
N=20
L=2.0
xmin = -1.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
    meshz[i] = xmin + 2*k*dx
end
umesh,region1 = chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,meshz)


--############################################### Execute meshing
--chiSurfaceMesherExecute();
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 21
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_graphite_pure.data")

src={}
for g=1,num_groups do
    src[g] = 0.0
end

chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics

phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,5, 5)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,20)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
if (master_export == nil) then
    --chiLBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
end
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)



fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

slice1 = chiFFInterpolationCreate(SLICE)
--chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.008,0.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0.0,0.0,1.0)
chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,0.0,-1.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,1.0,0.0,0.0)
chiFFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice1)
chiFFInterpolationExecute(slice1)

slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

chiFFInterpolationExportPython(slice1)
chiFFInterpolationExportPython(slice2)

if (chi_location_id == 0) then
    local handle = io.popen("python3 ZPFFI00.py")
    local handle = io.popen("python3 ZPFFI10.py")
    print("Execution completed")
end

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5e", maxval))

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[20])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end

if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")
end

