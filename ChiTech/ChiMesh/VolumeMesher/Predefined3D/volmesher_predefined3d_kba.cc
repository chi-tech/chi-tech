#include "volmesher_predefined3d.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell_polyhedron.h"


#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Gets the partition ID from a centroid for KBA-style partitioning.*/
int chi_mesh::VolumeMesherPredefined3D::
  GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid)
{
  auto handler = chi_mesh::GetCurrentHandler();

  int Px = handler->surface_mesher->partitioning_x;
  int Py = handler->surface_mesher->partitioning_y;


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
bool chi_mesh::VolumeMesherPredefined3D::
  IsRawCellNeighborToPartitionKBA(
  const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell)
{
  auto handler = chi_mesh::GetCurrentHandler();

  auto umesh = handler->unpartitionedmesh_stack.back();

  bool is_neighbor = false;
  for (const auto& face : lwcell.faces)
  {
    if (face.neighbor < 0) continue;
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
void chi_mesh::VolumeMesherPredefined3D::
  KBA(chi_mesh::UnpartitionedMesh* umesh,
      chi_mesh::MeshContinuumPtr grid)
{
  //======================================== Load up the vertices
  for (auto vert : umesh->vertices)
    grid->vertices.push_back(new chi_mesh::Vertex(*vert));

  chi_log.Log(LOG_0) << "Vertices loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  int loc_id = chi_mpi.location_id;

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

//    printf("[%d] Bla\n",loc_id);

    if (temp_cell->partition_id != chi_mpi.location_id)
    {
//      printf("[%d] Bla1\n",loc_id);
      if (IsRawCellNeighborToPartitionKBA(*raw_cell))
        grid->cells.push_back(temp_cell);
      else
        delete temp_cell;

//      printf("[%d] Bla2\n",loc_id);
    }
    else
    {
//      printf("[%d] Bla3\n",loc_id);
      auto polyh_cell = new chi_mesh::CellPolyhedron;
      polyh_cell->centroid = temp_cell->centroid;
      polyh_cell->global_id = temp_cell->global_id;
      polyh_cell->partition_id = temp_cell->partition_id;
      polyh_cell->material_id = temp_cell->material_id;

      polyh_cell->vertex_ids = raw_cell->vertex_ids;

      for (auto& raw_face : raw_cell->faces)
      {
        chi_mesh::CellFace newFace;

        if (raw_face.neighbor >= 0)
        {
          newFace.neighbor_id = raw_face.neighbor;
          newFace.has_neighbor = true;
        }

        newFace.vertex_ids = raw_face.vertex_ids;
        auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
        for (auto fvid : newFace.vertex_ids)
          vfc = vfc + *grid->vertices[fvid];
        newFace.centroid = vfc / newFace.vertex_ids.size();

        //Compute normal
//        auto va = *grid->vertices[newFace.vertex_ids[0]] - vfc;
//        auto vb = *grid->vertices[newFace.vertex_ids[1]] - vfc;
//
//        auto vn = va.Cross(vb);
//        newFace.normal = vn/vn.Norm();

        newFace.normal = chi_mesh::Normal(0.0,0.0,0.0);
        int last_vert_ind = newFace.vertex_ids.size()-1;
        for (int fv=0; fv<newFace.vertex_ids.size(); ++fv)
        {
          int fvid_m = newFace.vertex_ids[fv];
          int fvid_p = (fv == last_vert_ind)? newFace.vertex_ids[0] :
                       newFace.vertex_ids[fv+1];
          auto leg_m = *grid->vertices[fvid_m] - newFace.centroid;
          auto leg_p = *grid->vertices[fvid_p] - newFace.centroid;

          auto vn = leg_m.Cross(leg_p);

          newFace.normal = newFace.normal + vn.Normalized();
        }
        newFace.normal = (newFace.normal/newFace.vertex_ids.size()).Normalized();

        polyh_cell->faces.push_back(newFace);
      }

      grid->cells.push_back(polyh_cell);
    }//else
  }
}