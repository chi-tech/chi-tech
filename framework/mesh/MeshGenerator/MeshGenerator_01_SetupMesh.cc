#include "MeshGenerator.h"

#include "mesh/Cell/cell.h"
#include "graphs/GraphPartitioner.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

namespace chi_mesh
{
/**Executes the partitioner and configures the mesh as a real mesh.*/
std::shared_ptr<MeshContinuum>
MeshGenerator::SetupMesh(std::unique_ptr<UnpartitionedMesh> input_umesh_ptr)
{
  const auto& raw_cells = input_umesh_ptr->GetRawCells();
  const size_t num_raw_cells = raw_cells.size();

  ChiLogicalErrorIf(num_raw_cells == 0, "No cells in final input mesh");

  //============================================= Build cell graph and centroids
  typedef std::vector<uint64_t> CellGraphNode;
  typedef std::vector<CellGraphNode> CellGraph;
  CellGraph cell_graph;
  std::vector<chi_mesh::Vector3> cell_centroids;

  cell_graph.reserve(num_raw_cells);
  cell_centroids.reserve(num_raw_cells);
  for (const auto& raw_cell_ptr : raw_cells)
  {
    CellGraphNode cell_graph_node;
    for (auto& face : raw_cell_ptr->faces)
      if (face.has_neighbor) cell_graph_node.push_back(face.neighbor);

    cell_graph.push_back(cell_graph_node);
    cell_centroids.push_back(raw_cell_ptr->centroid);
  }

  //============================================= Execute partitioner
  auto cell_pids =
    partitioner_->Partition(cell_graph, cell_centroids, Chi::mpi.process_count);

  //============================================= Convert mesh
  auto grid_ptr = chi_mesh::MeshContinuum::New();

  grid_ptr->GetBoundaryIDMap() =
    input_umesh_ptr->GetMeshOptions().boundary_id_map;

  auto& vertex_subs = input_umesh_ptr->GetVertextCellSubscriptions();
  size_t cell_globl_id = 0;
  for (auto raw_cell : input_umesh_ptr->GetRawCells())
  {
    if (CellHasLocalScope(*raw_cell, cell_globl_id, vertex_subs, cell_pids))
    {
      auto cell = SetupCell(*raw_cell,
                            cell_globl_id,
                            cell_pids[cell_globl_id],
                            input_umesh_ptr->GetVertices());

      for (uint64_t vid : cell->vertex_ids_)
        grid_ptr->vertices.Insert(vid, input_umesh_ptr->GetVertices()[vid]);

      grid_ptr->cells.push_back(std::move(cell));
    }

    ++cell_globl_id;
  } // for raw_cell

  grid_ptr->SetGlobalVertexCount(input_umesh_ptr->GetVertices().size());

  //======================================== Concluding messages
  Chi::log.LogAllVerbose1()
    << "### LOCATION[" << Chi::mpi.location_id
    << "] amount of local cells=" << grid_ptr->local_cells.size();

  size_t total_local_cells = grid_ptr->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(&total_local_cells,
                &total_global_cells,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                Chi::mpi.comm);

  Chi::log.Log() << "MeshGenerator: Cells created = " << total_global_cells
                 << std::endl;

  return grid_ptr;
}

} // namespace chi_mesh