chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,"XMesh_Test3.obj",true)

region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);


chiSurfaceMesherCreate(SURFACEMESHER_TRIANGLE);
chiSurfaceMesherSetProperty(MAX_AREA,50/20/1.0)

chiSurfaceMesherExecute();

chiSurfaceMesherExportToObj("YMeshedSurf3.obj")

