chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,"Mesh/PlaneSurfaceMesh.obj",true)

region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);


chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);

chiSurfaceMesherExecute();

chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);
chiVolumeMesherSetProperty(EXTRUSION_LAYER,1.0,1,"Layer");

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);

chiSurfaceMesherSetProperty(PARTITION_X,CHI_PX)
chiSurfaceMesherSetProperty(PARTITION_Y,CHI_PY)
chiVolumeMesherSetProperty(PARTITION_Z,CHI_PZ)

chiVolumeMesherExecute()

surf_lv0 = chiSurfaceMeshCreate();

chiSurfaceMeshImportFromOBJFile(surf_lv0,"Mesh/LV_Cube.obj",false)

vol_lv0 = chiLogicalVolumeCreate(SURFACE,surf_lv0)

chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol_lv0,0)

materials = {}
materials[1] = chiPhysicsAddMaterial("Material 1")
