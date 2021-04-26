#include "volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"


#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Gets the partition ID from a centroid for KBA-style partitioning.*/
int chi_mesh::VolumeMesherPredefinedUnpartitioned::
  GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid)
{
  auto handler = chi_mesh::GetCurrentHandler();

  int Px = handler->volume_mesher->options.partition_x;
  int Py = handler->volume_mesher->options.partition_y;

  chi_mesh::Cell temp_cell(chi_mesh::CellType::GHOST);
  temp_cell.centroid = centroid;

  auto xyz = GetCellXYZPartitionID(&temp_cell);

  int nxi = std::get<0>(xyz);
  int nyi = std::get<1>(xyz);
  int nzi = std::get<2>(xyz);

  return nzi*Px*Py + nyi*Px + nxi;
}

//###################################################################
/**Determines if a chi_mesh::UnpartitionedMesh::LightWeightCell is a
 * neighbor to the current partition for KBA-style partitioning.
 * This method loops over the faces of the lightweight cell and
 * determines the partition-id of each the neighbors. If the neighbor
 * has a partition id equal to that of the current process then
 * it means this reference cell is a neighbor.*/
bool chi_mesh::VolumeMesherPredefinedUnpartitioned::
IsRawCellNeighborToPartitionKBA(
  const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell)
{
  auto handler = chi_mesh::GetCurrentHandler();

  auto umesh = handler->unpartitionedmesh_stack.back();

  bool is_neighbor = false;
  for (const auto& face : lwcell.faces)
  {
    if (not face.has_neighbor) continue;
    auto adj_cell = umesh->raw_cells[face.neighbor];
    int partition_id = GetPartitionIDFromCentroid(adj_cell->centroid);
    if (partition_id == chi_mpi.location_id)
    {
      is_neighbor = true;
      break;
    }
  }

  return is_neighbor;
}

//###################################################################
/** Applies KBA-style partitioning to the mesh.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
KBA(chi_mesh::UnpartitionedMesh* umesh,
    chi_mesh::MeshContinuumPtr& grid)
{
  //======================================== Load up the vertices
  for (auto& vert : umesh->vertices)
    grid->vertices.push_back(vert);

  chi_log.Log(LOG_0) << "Vertices loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  //======================================== Load up the cells
  int global_id=-1;
  for (auto raw_cell : umesh->raw_cells)
  {
    ++global_id;
    auto temp_cell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
    temp_cell->centroid = raw_cell->centroid;
    temp_cell->global_id = global_id;
    temp_cell->partition_id = GetPartitionIDFromCentroid(temp_cell->centroid);
    temp_cell->material_id = raw_cell->material_id;

    if (temp_cell->partition_id != chi_mpi.location_id)
    {
      if (IsRawCellNeighborToPartitionKBA(*raw_cell))
        grid->cells.push_back(temp_cell);
      else
        delete temp_cell;
    }
    else
    {
      if (raw_cell->type == chi_mesh::CellType::SLAB)
        AddSlabToGrid(*raw_cell,*temp_cell,*grid);
      else if (raw_cell->type == chi_mesh::CellType::POLYGON)
        AddPolygonToGrid(*raw_cell,*temp_cell,*grid);
      else if (raw_cell->type == chi_mesh::CellType::POLYHEDRON)
        AddPolyhedronToGrid(*raw_cell,*temp_cell,*grid);
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "Unsupported cell encountered in "
             "chi_mesh::VolumeMesherPredefinedUnpartitioned";
        exit(EXIT_FAILURE);
      }
    }//else
  }
}