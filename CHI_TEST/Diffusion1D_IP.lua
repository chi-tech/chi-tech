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
line_mesh = chiLineMeshCreateFromArray(mesh)


region1 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);


--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

chiVolumeMesherSetProperty(PARTITION_Z,2)


--############################################### Execute meshing
chiSurfaceMesherExecute();
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
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLD_MIP);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-6)



--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
fftemp,count = chiGetFieldFunctionList(phys1)
--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,0,DIFFUSION_DIRICHLET,0.0)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,1,DIFFUSION_DIRICHLET,0.0)

chiDiffusionExecute(phys1)

line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,0.1,0.0,-1.0)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT,0.1,0.0, 1.0)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(line0)
chiFFInterpolationExecute(line0)
--chiFFInterpolationExportPython(line0)

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.10f", maxval))
--
if ((chi_location_id == 0) and (master_export == nil)) then
    chiFFInterpolationExportPython(line0)
    local handle = io.popen("python ZLFFI00.py")
end

--chiExportFieldFunctionToVTK(fftemp,"ZPhi1D")