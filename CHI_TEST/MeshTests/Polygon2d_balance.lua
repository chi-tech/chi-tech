if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

chiDecomposeSurfaceMeshPxPy(newSurfMesh,4,4)

chiSurfaceMeshExtractOpenEdgesToObj(newSurfMesh,"ZEdgesTest.obj")