#include "volmesher_extruder.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/** Creates actual z-levels for the input layer specification.*/
void chi_mesh::VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  //================================================== Create default layers if no
  //                                                   input layers are provided
  if (input_layers.empty())
  {
    chi_log.Log(LOG_0WARNING)
      << "VolumeMesherExtruder: No extrusion layers have been specified. "
      << "A default single layer will be used with height 1.0 and a single "
      << "subdivision.";
    double dz = 1.0/default_layer_count;
    for (int k=0; k<=default_layer_count; k++)
    {
      vertex_layers.push_back(k*dz);
    }
  }
  else
  {
    double last_z=0.0;
    vertex_layers.push_back(last_z);

    for (const auto& input_layer : input_layers)
    {
      double dz = input_layer.height/input_layer.sub_divisions;

      for (int k=0; k<input_layer.sub_divisions; k++)
      {
        last_z += dz;
        vertex_layers.push_back(last_z);
      }
    }
  }

  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Total number of cell layers is "
    << vertex_layers.size()-1;
}

//###################################################################
/**Projects a centroid to an extruded equivalent layer.*/
chi_mesh::Vector3 chi_mesh::VolumeMesherExtruder::
  ProjectCentroidToLevel(const chi_mesh::Vector3 &centroid,
                         const size_t level)
{
  double z_location = 0.5 * (vertex_layers[level] + vertex_layers[level+1]);

  auto centroid_projected = centroid;
  centroid_projected.z = z_location;

  return centroid_projected;
}

//###################################################################
/**Computes a cell's partition id based on a centroid.*/
int chi_mesh::VolumeMesherExtruder::
  GetCellKBAPartitionIDFromCentroid(chi_mesh::Vector3& centroid)
{
  int px = options.partition_x;
  int py = options.partition_y;

  chi_mesh::Cell n_gcell(CellType::GHOST, CellType::GHOST);
  n_gcell.centroid = centroid;

  auto xyz_partition_indices = GetCellXYZPartitionID(&n_gcell);

  int nxi = std::get<0>(xyz_partition_indices);
  int nyi = std::get<1>(xyz_partition_indices);
  int nzi = std::get<2>(xyz_partition_indices);

  return nzi*px*py + nyi*px + nxi;
}


//###################################################################
/**Determines if a template cell is in the current partition or a
 * direct neighbor to the current partition.*/
bool chi_mesh::VolumeMesherExtruder::
  HasLocalScope(
    const chi_mesh::Cell& template_cell,
    const chi_mesh::MeshContinuum& template_continuum,
    size_t z_level)
{
  //======================================== Check if the template cell
  //                                         is in the current partition
  {
    auto& centroid = template_cell.centroid;
    auto projected_centroid = ProjectCentroidToLevel(centroid, z_level);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == chi_mpi.location_id)
      return true;
  }

  const size_t last_z_level = vertex_layers.size() - 1;
  const size_t z_level_below = z_level - 1;
  const size_t z_level_above = z_level + 1;

  //======================================== Build z-levels to search
  std::vector<size_t> z_levels_to_search;
  z_levels_to_search.reserve(3);

                                   z_levels_to_search.push_back(z_level);
  if (z_level != 0)                z_levels_to_search.push_back(z_level_below);
  if (z_level != (last_z_level-1)) z_levels_to_search.push_back(z_level_above);

  //======================================== Search template cell's
  //                                         lateral neighbors
  const auto& vertex_subs =
    template_unpartitioned_mesh->vertex_cell_subscriptions;
  for (uint64_t vid : template_cell.vertex_ids)
    for (uint64_t cid : vertex_subs[vid])
    {
      if (cid == template_cell.local_id) continue;

      auto& candidate_cell =
        template_continuum.local_cells[cid];
      auto& cc_centroid = candidate_cell.centroid;

      for (size_t z : z_levels_to_search)
      {
        auto projected_centroid = ProjectCentroidToLevel(cc_centroid, z);
        int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
        if (pid == chi_mpi.location_id)
          return true;
      }
    }//for cid

  //======================================== Search template cell's
  //                                         longitudinal neighbors
  for (size_t z : z_levels_to_search)
  {
    if (z == z_level) continue;

    auto projected_centroid = ProjectCentroidToLevel(template_cell.centroid, z);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == chi_mpi.location_id)
      return true;
  }//for z

  return false;
}

//###################################################################
/**Makes an extruded cell from a template cell.*/
chi_mesh::Cell* chi_mesh::VolumeMesherExtruder::
  MakeExtrudedCell(const chi_mesh::Cell &template_cell,
                   const chi_mesh::MeshContinuum &grid,
                   size_t z_level,
                   uint64_t cell_global_id,
                   int partition_id,
                   size_t num_template_cells)
{
  const size_t tc_num_verts = template_cell.vertex_ids.size();

  //========================================= Create polyhedron
  auto cell = new chi_mesh::Cell(CellType::POLYHEDRON, CellType::POLYHEDRON);
  cell->global_id    = cell_global_id;
  //cell->local_id set when added to mesh
  cell->partition_id = partition_id;
  cell->centroid = ProjectCentroidToLevel(template_cell.centroid, z_level);

  //========================================= Populate cell v-indices
  cell->vertex_ids.reserve(2*tc_num_verts);
  for (auto tc_vid : template_cell.vertex_ids)
    cell->vertex_ids.push_back(tc_vid + z_level*node_z_index_incr);

  for (auto tc_vid : template_cell.vertex_ids)
    cell->vertex_ids.push_back(tc_vid + (z_level+1)*node_z_index_incr);


  //========================================= Create side faces
  for (auto& face : template_cell.faces)
  {
    chi_mesh::CellFace newFace;

    newFace.vertex_ids.resize(4,-1);
    newFace.vertex_ids[0] = face.vertex_ids[0] + z_level*node_z_index_incr;
    newFace.vertex_ids[1] = face.vertex_ids[1] + z_level*node_z_index_incr;
    newFace.vertex_ids[2] = face.vertex_ids[1] + (z_level+1)*node_z_index_incr;
    newFace.vertex_ids[3] = face.vertex_ids[0] + (z_level+1)*node_z_index_incr;

    //Compute centroid
    const auto& v0 = grid.vertices[newFace.vertex_ids[0]];
    const auto& v1 = grid.vertices[newFace.vertex_ids[1]];
    const auto& v2 = grid.vertices[newFace.vertex_ids[2]];
    const auto& v3 = grid.vertices[newFace.vertex_ids[3]];

    chi_mesh::Vertex vfc = (v0+v1+v2+v3)/4.0;
    newFace.centroid = vfc;

    //Compute normal
    chi_mesh::Vector3 va = v0 - vfc;
    chi_mesh::Vector3 vb = v1 - vfc;

    chi_mesh::Vector3 vn = va.Cross(vb);

    newFace.normal = (vn/vn.Norm());

    //Set neighbor
    //The side connections have the same connections as the
    //template cell + the z_level specifiers of the layer.
    if (face.has_neighbor)
    {
      newFace.neighbor_id = face.neighbor_id +
                            z_level*num_template_cells;
      newFace.has_neighbor = true;
    }
    else
      newFace.neighbor_id = face.neighbor_id;

    cell->faces.push_back(newFace);
  } //for side faces

  //========================================= Create top and bottom faces
  //=============================== Bottom face
  {
    chi_mesh::CellFace newFace;
    //Vertices
    auto vfc = chi_mesh::Vertex(0.0,0.0,0.0);
    newFace.vertex_ids.reserve(tc_num_verts);
    for (int tv=(static_cast<int>(tc_num_verts)-1); tv>=0; tv--)
    {
      newFace.vertex_ids.push_back(template_cell.vertex_ids[tv]
                                   + z_level*node_z_index_incr);
      const auto& v = grid.vertices[newFace.vertex_ids.back()];
      vfc = vfc + v;
    }

    //Compute centroid
    vfc = vfc/static_cast<double>(tc_num_verts);
    newFace.centroid = vfc;

    //Compute normal
    auto va = grid.vertices[newFace.vertex_ids[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal = vn/vn.Norm();

    //Set neighbor
    if (z_level==0)
    {
      newFace.neighbor_id = 0;
      newFace.has_neighbor = false;
    }
    else
    {
      newFace.neighbor_id = template_cell.local_id +
        (z_level-1)*num_template_cells;
      newFace.has_neighbor = true;
    }

    cell->faces.push_back(newFace);
  }

  //=============================== Top face
  {
    chi_mesh::CellFace newFace;
    //Vertices
    auto vfc = chi_mesh::Vertex(0.0,0.0,0.0);
    newFace.vertex_ids.reserve(tc_num_verts);
    for (auto tc_vid : template_cell.vertex_ids)
    {
      newFace.vertex_ids.push_back(tc_vid + (z_level+1)*node_z_index_incr);
      const auto& v = grid.vertices[newFace.vertex_ids.back()];
      vfc = vfc + v;
    }

    //Compute centroid
    vfc = vfc/static_cast<double>(tc_num_verts);
    newFace.centroid = vfc;

    //Compute normal
    auto va = grid.vertices[newFace.vertex_ids[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal = vn/vn.Norm();

    //Set neighbor
    if (z_level==(vertex_layers.size()-2))
    {
      newFace.neighbor_id = 0;
      newFace.has_neighbor = false;
    }
    else
    {
      newFace.neighbor_id = template_cell.local_id +
        (z_level+1)*num_template_cells;
      newFace.has_neighbor = true;
    }

    cell->faces.push_back(newFace);
  }

  return cell;
}