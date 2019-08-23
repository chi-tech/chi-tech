if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);

NZ=10
chiVolumeMesherSetProperty(EXTRUSION_LAYER,1.0,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--chiVolumeMesherSetProperty(EXTRUSION_LAYER,0.2,NZ,"Charlie");
--
chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,-0.5)
chiSurfaceMesherSetProperty(CUT_Y,0.0)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-0.5,0.5)
vol2 = chiLogicalVolumeCreate(RPP,-0.25,0.25,-0.25,0.25,-0.25,0.25)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol2,2)

chiPhysicsAddMaterial("Air");
chiPhysicsAddMaterial("Graphite");
chiPhysicsAddMaterial("HDPE");

chiRegionExportMeshToVTK(region1,"ZHMesh")