chiMeshHandlerCreate()

--chiUnpartitionedMeshFromVTU("CHI_TEST/UnpartTest/Z_VTK_Mesh_0.vtu")
chiUnpartitionedMeshFromEnsightGold("CHI_TEST/UnpartTest/Sphere.case")

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED3D)

chiVolumeMesherSetProperty(CUTS_X,0.0)
chiVolumeMesherSetProperty(CUTS_Y,0.0)

chiVolumeMesherSetProperty(VOLUMEPARTITION_X,2)
chiVolumeMesherSetProperty(VOLUMEPARTITION_Y,2)




chiSurfaceMesherExecute()
chiVolumeMesherExecute()

chiRegionExportMeshToVTK(region1,"Mesh")


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[1],SCALAR_VALUE,SINGLE_VALUE,1.0)


chiPhysicsMaterialAddProperty(materials[2],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[2],SCALAR_VALUE,SINGLE_VALUE,1.0)



--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLD_MIP);
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

chiExportFieldFunctionToVTK(fftemp[1],"ZPhi")