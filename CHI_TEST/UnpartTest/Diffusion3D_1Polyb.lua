if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



---############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
meshz={}
N=40
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


--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)



--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-6)

--############################################### Set boundary conditions
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIRICHLET,0.1)

--############################################### Set boundary conditions
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,5,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,6,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,0,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,2,DIFFUSION_VACUUM)
--chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,3,DIFFUSION_VACUUM)

--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
fftemp,count = chiGetFieldFunctionList(phys1)


chiDiffusionExecute(phys1)

--############################################### Set derived geometry
slice1 = chiFFInterpolationCreate(SLICE)
--chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice1,SLICE_POINT,0.008,0.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_BINORM,0.0,0.0,1.0)
chiFFInterpolationSetProperty(slice1,SLICE_TANGENT,0.0,-1.0,0.0)
chiFFInterpolationSetProperty(slice1,SLICE_NORMAL,1.0,0.0,0.0)
chiFFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fftemp[1])


slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp[1])

line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.0,0.025)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.0,0.025)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(slice1)
chiFFInterpolationExecute(slice1)


chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)


chiFFInterpolationInitialize(line0)
chiFFInterpolationExecute(line0)


ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

if (master_export == nil) then
    chiFFInterpolationExportPython(slice1)
    chiFFInterpolationExportPython(slice2)
    chiFFInterpolationExportPython(line0)
end

if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python3 ZPFFI00.py")
    local handle = io.popen("python3 ZPFFI10.py")
    local handle = io.popen("python3 ZLFFI20.py")
    print("Execution completed")
end

if (master_export == nil) then
    chiExportFieldFunctionToVTK(fftemp,"ZPhi")
end
