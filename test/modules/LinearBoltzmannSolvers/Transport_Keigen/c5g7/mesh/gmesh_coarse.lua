--############################################### Setup mesh
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromMshFormat("mesh/2D_c5g7_coarse.msh",true)

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED,umesh);

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();
