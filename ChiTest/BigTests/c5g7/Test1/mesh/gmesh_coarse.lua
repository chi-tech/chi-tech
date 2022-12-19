--############################################### Setup mesh
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromMshFormat("mesh/2D_c5g7_coarse.msh",true)

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED,umesh);

--chiVolumeMesherSetKBAPartitioningPxPyPz(2,2,1)
--chiVolumeMesherSetKBACutsX({27.72})
--chiVolumeMesherSetKBACutsY({27.72})
--chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)


--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();
