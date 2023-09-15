/**\defgroup MeshGenerator Mesh Generators
*
We split a `MeshGenerator`'s execution into a
phase that generates an unpartitioned mesh and a phase that then converts
this mesh into partitioned `chi_mesh::MeshContinuum` (with both steps
customizable). The phase that creates the `MeshContinuum` object can be hooked
up to a partitioner that can also be designed to be pluggable.

## Example A
\code
nodes = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0}
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
chi_mesh.MeshGenerator.Execute(meshgen1)

chiMeshHandlerExportMeshToVTK("ZMeshTest")
\endcode

In this example we created a set of nodes (monotonically increasing in value)
for use with the
\ref chi_mesh__OrthogonalMeshGenerator. We supplied the same set twice meaning the generator
will build a 2D mesh.

We did not specify a partitioner and therefore the generator will use the
\ref chi__PETScGraphPartitioner with `type="parmetis"` by default. Using 8
processes, we can see the mesh and it's partitioning below.

\image html framework/chi_mesh/MeshGenerators/ExampleA.png width=500px

## Example B
\code
meshgen1 = chi_mesh.FromFileMeshGenerator.Create({ filename="TriangleMesh2x2.obj" })
chi_mesh.MeshGenerator.Execute(meshgen1)

chiMeshHandlerExportMeshToVTK("ZMeshTest")
\endcode

In this example we created a mesh by reading it from a file, using the
\ref chi_mesh__FromFileMeshGenerator. The types of file-types we can support is
ever growing. At the time of writing this we support the following formats:
- `.obj` Wavefront
- `.msh` gmesh,
- `.e` ExodusII,
- `.vtu` VTK Unstructured grid,
- `.pvtu` Pieced VTK Unstructured grid,
- `.case` Ensight Gold

Using 8 processes, we can see the mesh and it's partitioning below.

\image html framework/chi_mesh/MeshGenerators/ExampleB.png width=500px

## Example C
\code
meshgen1 = chi_mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create({ filename=TriangleMesh2x2.obj" }),
  },
  layers = {{z=1.1, n=2},    -- First layer - 2 sub-layers
            {z=2.1, n=3}},   -- Second layer - 3 sub-layers
})
chi_mesh.MeshGenerator.Execute(meshgen1)

chiMeshHandlerExportMeshToVTK("ZMeshTest")
\endcode

In this example we use an \ref chi_mesh__ExtruderMeshGenerator with another
mesh generator as an input to create an extruded mesh. Using 8 processes, we can
see the mesh and it's partitioning below.

\image html framework/chi_mesh/MeshGenerators/ExampleC.png width=500px

## Using different Partitioners
ChiTech now has a set of Graph Partitioners (i.e. based off `GraphPartitioner`)
that support different forms of partitioning, for example we have:
- \ref chi__LinearGraphPartitioner, a very basic partitioner used for
during the preparation of simulation meshes.
- \ref chi__PETScGraphPartitioner, a flexible partitioner that can use all the
partitioner options available in PETSc (defaults to using `"parmetis"`).
- \ref chi__KBAGraphPartitioner, the classical Neutron Transport KBA parallel
partitioning with an overlaid orthogonal layout.

An example of changing the partitioning to PETSc's `"average"` option is shown
below:
\code
meshgen1 = chi_mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename="resources/TestMeshes/TriangleMesh2x2.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}},
  partitioner = chi.PETScGraphPartitioner.Create({type="average"})
})
chi_mesh.MeshGenerator.Execute(meshgen1)

chiMeshHandlerExportMeshToVTK("ZMeshTest")
\endcode

\image html framework/chi_mesh/MeshGenerators/ParExample1.png width=500px

Another example using the `KBAGraphPartitioner` is shown below
\code
meshgen1 = chi_mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename="resources/TestMeshes/TriangleMesh2x2.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}},
  partitioner = chi.KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=2,
    xcuts = {0.0}, ycuts = {0.0}, zcuts = {1.1}
  })
})
chi_mesh.MeshGenerator.Execute(meshgen1)

chiMeshHandlerExportMeshToVTK("ZMeshTest")
\endcode

\image html framework/chi_mesh/MeshGenerators/ParExample2.png width=500px

\ingroup LuaMesh*/
