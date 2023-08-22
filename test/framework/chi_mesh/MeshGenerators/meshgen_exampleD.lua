meshgen1 = chi_mesh.MeshGenerator.Create
({
  inputs =
  {
    [0] = chi_mesh.FromFileMeshGenerator.Create
          ({
            filename="TriangleMesh2x2.obj"
          }),
    [1] = chi_mesh.ExtruderMeshGenerator.Create
          ({
            layers = {{z=1.1, n=2}, {z=2.1, n=3}}
          })
  }
})
chi_mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
