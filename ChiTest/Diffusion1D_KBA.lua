print("############################################### LuaTest")
--dofile(CHI_LIBRARY)

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=100
L=2.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)

print(PARTITION_TYPE,KBA_STYLE_XYZ,CUT_Z,PARTITION_Z)
chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)
-- chiVolumeMesherSetProperty(CUTS_Z,L/2.0)
chiVolumeMesherSetProperty(PARTITION_Z,2)


chiVolumeMesherExecute();


--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
chiVolumeMesherSetupOrthogonalBoundaries()


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

chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,OrthoBoundaryID.ZMIN,DIFFUSION_VACUUM)
chiDiffusionSetProperty(phys1,BOUNDARY_TYPE,OrthoBoundaryID.ZMAX,DIFFUSION_VACUUM)


--############################################### Initialize Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

fftemp,count = chiGetFieldFunctionList(phys1)
ffi0 = chiFFInterpolationCreate(LINE)
curffi = ffi0;
chiFFInterpolationSetProperty(curffi,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(curffi,LINE_SECONDPOINT,0.0,0.0, 2.0+xmin)
chiFFInterpolationSetProperty(curffi,LINE_NUMBEROFPOINTS, 1000)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

if (chi_location_id == 0 and master_export == nil) then
    chiFFInterpolationExportPython(ffi0)

    local handle = io.popen("python ZLFFI00.py")
end
