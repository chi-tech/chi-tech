-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.51187 and 1.42458e-03
num_procs = 4
--Unstructured mesh




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromWavefrontOBJ("tests/MeshIO/MeshWavefront3.obj")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh);

chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Add materials
materials = {}
materials["0"] = chiPhysicsAddMaterial("Test Material");
materials["1"] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials["0"],SCALAR_VALUE, "D")
chiPhysicsMaterialAddProperty(materials["0"],SCALAR_VALUE, "q")
chiPhysicsMaterialSetProperty(materials["0"],"D",SINGLE_VALUE,1.0)
chiPhysicsMaterialSetProperty(materials["0"],"q",SINGLE_VALUE,100.0)

chiPhysicsMaterialAddProperty(materials["1"],SCALAR_VALUE, "D")
chiPhysicsMaterialAddProperty(materials["1"],SCALAR_VALUE, "q")
chiPhysicsMaterialSetProperty(materials["1"],"D",SINGLE_VALUE,10.0)
chiPhysicsMaterialSetProperty(materials["1"],"q",SINGLE_VALUE,0.0)

--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLC")
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-8)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,"boundary_type","BndryTop","dirichlet",0.0)
chiDiffusionSetProperty(phys1,"boundary_type","BndryBottom","dirichlet",0.0)
chiDiffusionSetProperty(phys1,"boundary_type","BndryLeftA","dirichlet",0.0)
chiDiffusionSetProperty(phys1,"boundary_type","BndryLeftB","dirichlet",1.0)
chiDiffusionSetProperty(phys1,"boundary_type","BndryRight","dirichlet",0.0)

--############################################### Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = chiSolverGetFieldFunctionList(phys1)

--############################################### Exports
if master_export == nil then
    chiMeshHandlerExportMeshToVTK("ZObjMesh2")
    chiExportMultiFieldFunctionToVTK(fftemp,"ZPhiObjMesh2")
end
