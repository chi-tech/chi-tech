/**\page DevManMeshes The structure of meshes

\section devman4_sec0 General structure of meshes

Meshes and mesh operations are all handled under the umbrella of a
chi_mesh::MeshHandler. If any mesh operation is performed it should push
newly created objects to a stack in a handler. The figure below shows
a generalized architecture.

\image html MeshOverview.png "Overview of the mesh hierarchy" width=700px

The final level of a chi_mesh::Cell object is a chi_mesh::MeshContinuum object.
On this object it will be placed inside a global stack, meaning all processes
will have the same stack (albeit some locations might have empty references for
cells not locally owned). Partitioning of meshes is facilitated by storing
the full data structure of a cell only if the cell is locally owned, and by
populating the local_cell_glob_index array of indices mapped from local ids
to global ids. This vector stores the global indices of cells that are locally
owned. See Figure below:

\image html CellMapping.png "Cell mapping logic" width=700px
 */