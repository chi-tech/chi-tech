chiMeshHandlerCreate()

chiUnpartitionedMeshFromEnsightGold("/Users/janv4/Desktop/GoogleDrive/Temp/DMDMeshes/Mesh3c.case",2.0)

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED3D)

chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)



chiSurfaceMesherExecute()
chiVolumeMesherExecute()

chiRegionExportMeshToVTK(region1,"Mesh")


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Air");
materials[2] = chiPhysicsAddMaterial("Ground");
materials[3] = chiPhysicsAddMaterial("ROI");
materials[4] = chiPhysicsAddMaterial("Shield");
materials[5] = chiPhysicsAddMaterial("Source");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[3],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[4],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[5],TRANSPORT_XSECTIONS)

--chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
--chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
--chiPhysicsMaterialAddProperty(materials[3],ISOTROPIC_MG_SOURCE)
--chiPhysicsMaterialAddProperty(materials[4],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[5],ISOTROPIC_MG_SOURCE)


num_groups = 2
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_air50RH.data")
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_air50RH.data")
--chiPhysicsMaterialSetProperty(materials[3],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_air50RH.data")
--chiPhysicsMaterialSetProperty(materials[4],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_air50RH.data")
--chiPhysicsMaterialSetProperty(materials[5],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_air50RH.data")

for i=1,5 do
    chiPhysicsMaterialSetProperty(materials[i],TRANSPORT_XSECTIONS,
            SIMPLEXS1,num_groups,1.0,0.1)
end


src={}
for g=1,num_groups do
    src[g] = 0.0
end

--chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1]=1.0e10
chiPhysicsMaterialSetProperty(materials[5],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics

phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad  = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,6, 6)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
if (master_export == nil) then
    --chiLBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
end
--chiLBSGroupsetSetMaxIterations(phys1,cur_gs,10)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0e12/4.0/math.pi;
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)



fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
--slices = {}
--for k=1,count do
--    slices[k] = chiFFInterpolationCreate(SLICE)
--    chiFFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,0.8001)
--    chiFFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
--    --chiFFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
--    --chiFFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
--    --chiFFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
--    chiFFInterpolationInitialize(slices[k])
--    chiFFInterpolationExecute(slices[k])
--    chiFFInterpolationExportPython(slices[k])
--end

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
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
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[2])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

slice1 = chiFFInterpolationCreate(SLICE)

--chiFFInterpolationSetProperty(slice1,SLICE_POINT,10000.0,10000.0,79900.0)
--chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0,0.0,-1.0)
--chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,100.0,135.71,0.0)
--chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,135.71,-100.0,0.0)

chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0,0.0,-1.0)
chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,1,0,0)
chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,0.0,-1,0.0)

chiFFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice1)
chiFFInterpolationExecute(slice1)

if (master_export == nil) then
    chiFFInterpolationExportPython(slice1)
end

if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI20.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end

if (master_export == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhi3D","Phi")
end


