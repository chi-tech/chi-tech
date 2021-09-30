print("############################################### LuaTest")
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=100
L=2.0
xmin = -1.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
umesh,region1 = chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)


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
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-4)

--############################################### Set boundary conditions
--chiDiffusionSetProperty(phys1,"boundary_type",0,"reflecting",1.0)
--chiDiffusionSetProperty(phys1,"boundary_type",1,"vacuum",2.0)
--chiDiffusionSetProperty(phys1,"boundary_type",2,"reflecting",3.0)
--chiDiffusionSetProperty(phys1,"boundary_type",3,"vacuum",4.0)

--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
fftemp,count = chiGetFieldFunctionList(phys1)
chiDiffusionExecute(phys1)

slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)


line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.1,0.0)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.1,0.0)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

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
    chiFFInterpolationExportPython(slice2)
    chiFFInterpolationExportPython(line0)
    local handle = io.popen("python3 ZPFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")

    chiExportFieldFunctionToVTK(fftemp,"ZPhi")
end
