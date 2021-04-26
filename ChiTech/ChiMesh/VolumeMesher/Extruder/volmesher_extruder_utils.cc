#include "volmesher_extruder.h"

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Computes the partition id of a template cell's projection onto 3D.*/
chi_mesh::Vector3
chi_mesh::VolumeMesherExtruder::ComputeTemplateCell3DCentroid(
  const chi_mesh::CellPolygon& n_template_cell,
  const chi_mesh::MeshContinuum& template_continuum,
  int z_level_begin,int z_level_end)
{
  chi_mesh::Vector3 n_centroid_precompd;
  size_t tc_num_verts = n_template_cell.vertex_ids.size();

  for (auto tc_vid : n_template_cell.vertex_ids)
  {
    auto temp_vert = template_continuum.vertices[tc_vid];
    temp_vert.z = vertex_layers[z_level_begin];
    n_centroid_precompd = n_centroid_precompd + temp_vert;
  }

  for (auto tc_vid : n_template_cell.vertex_ids)
  {
    auto temp_vert = template_continuum.vertices[tc_vid];
    temp_vert.z = vertex_layers[z_level_end];
    n_centroid_precompd = n_centroid_precompd + temp_vert;
  }
  n_centroid_precompd = n_centroid_precompd/(2*(double)tc_num_verts);

  return n_centroid_precompd;
}

//###################################################################
/**Computes a cell's partition id based on a centroid.*/
int chi_mesh::VolumeMesherExtruder::
GetCellPartitionIDFromCentroid(chi_mesh::Vector3& centroid)
{
  int px = options.partition_x;
  int py = options.partition_y;

  chi_mesh::Cell n_gcell(chi_mesh::CellType::GHOST);
  n_gcell.centroid = centroid;

  auto xyz_partition_indices = GetCellXYZPartitionID(&n_gcell);

  int nxi = std::get<0>(xyz_partition_indices);
  int nyi = std::get<1>(xyz_partition_indices);
  int nzi = std::get<2>(xyz_partition_indices);

  return nzi*px*py + nyi*px + nxi;
}

//###################################################################
/**Determines if a template cell is neighbor to the current partition.*/
bool chi_mesh::VolumeMesherExtruder::
IsTemplateCellNeighborToThisPartition(
  const chi_mesh::CellPolygon& template_cell,
  const chi_mesh::MeshContinuum& template_continuum,
  int z_level,int tc_index)
{
  int iz = z_level;
  int tc = tc_index;

  //========================= Loop over template cell neighbors
  //                          for side neighbors
  bool is_neighbor_to_partition = false;
  for (auto& tc_face : template_cell.faces)
  {
    if (tc_face.has_neighbor)
    {
      auto n_template_cell = (chi_mesh::CellPolygon&)(
        template_continuum.local_cells[tc_face.neighbor_id]);

      auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
        n_template_cell, template_continuum, iz, iz+1);

      int n_gcell_partition_id =
        GetCellPartitionIDFromCentroid(n_centroid_precompd);

      if (n_gcell_partition_id == chi_mpi.location_id)
      {
        is_neighbor_to_partition = true;
        break;
      }
    }//if neighbor not border
  }//for neighbors

  //========================= Now look at bottom neighbor
  //Bottom face
  if (iz != 0)
  {
    auto n_template_cell = (chi_mesh::CellPolygon&)
      (template_continuum.local_cells[tc]);

    auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
      n_template_cell, template_continuum, iz-1, iz);

    int n_gcell_partition_id =
      GetCellPartitionIDFromCentroid(n_centroid_precompd);

    if (n_gcell_partition_id == chi_mpi.location_id)
      is_neighbor_to_partition = true;
  }//if neighbor not border

  //========================= Now look at top neighbor
  //Top Face
  if (iz != (vertex_layers.size()-2))
  {
    auto n_template_cell = (chi_mesh::CellPolygon&)
      (template_continuum.local_cells[tc]);

    auto n_centroid_precompd = ComputeTemplateCell3DCentroid(
      n_template_cell, template_continuum, iz+1, iz+2);

    int n_gcell_partition_id =
      GetCellPartitionIDFromCentroid(n_centroid_precompd);

    if (n_gcell_partition_id == chi_mpi.location_id)
      is_neighbor_to_partition = true;
  }//if neighbor not border

  return is_neighbor_to_partition;
}