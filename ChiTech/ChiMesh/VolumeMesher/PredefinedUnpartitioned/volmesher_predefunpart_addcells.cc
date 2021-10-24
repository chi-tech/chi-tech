#include "volmesher_predefunpart.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Adds a slab to the grid from a light-weight cell.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  AddSlabToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell &raw_cell,
    const chi_mesh::Cell &temp_cell,
    chi_mesh::MeshContinuum &grid)
{
  auto slab_cell = new chi_mesh::Cell(CellType::SLAB, CellType::SLAB);
  slab_cell->centroid = temp_cell.centroid;
  slab_cell->global_id = temp_cell.global_id;
  slab_cell->partition_id = temp_cell.partition_id;
  slab_cell->material_id = temp_cell.material_id;

  slab_cell->vertex_ids = raw_cell.vertex_ids;

  int fc=-1;
  for (auto& raw_face : raw_cell.faces)
  {
    ++fc;
    chi_mesh::CellFace newFace;

    newFace.has_neighbor = raw_face.has_neighbor;
    newFace.neighbor_id = raw_face.neighbor;

    newFace.vertex_ids = raw_face.vertex_ids;
    auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids)
      vfc = vfc + grid.vertices[fvid];
    newFace.centroid = vfc / double(newFace.vertex_ids.size());

    // A slab face is very easy. If it is the first face
    // the normal is -khat. If it is the second face then
    // it is +khat.
    if (fc == 0)
      newFace.normal = chi_mesh::Vector3(0.0,0.0,-1.0);
    else
      newFace.normal = chi_mesh::Vector3(0.0,0.0, 1.0);

    slab_cell->faces.push_back(newFace);
  }

  grid.cells.push_back(slab_cell);
}

//###################################################################
/**Adds a polygon to the grid from a light-weight cell.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  AddPolygonToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell &raw_cell,
    const chi_mesh::Cell &temp_cell,
    chi_mesh::MeshContinuum &grid)
{
  auto poly_cell = new chi_mesh::Cell(CellType::POLYGON, raw_cell.sub_type);
  poly_cell->centroid = temp_cell.centroid;
  poly_cell->global_id = temp_cell.global_id;
  poly_cell->partition_id = temp_cell.partition_id;
  poly_cell->material_id = temp_cell.material_id;

  poly_cell->vertex_ids = raw_cell.vertex_ids;

  for (auto& raw_face : raw_cell.faces)
  {
    chi_mesh::CellFace newFace;

    newFace.has_neighbor = raw_face.has_neighbor;
    newFace.neighbor_id = raw_face.neighbor;

    newFace.vertex_ids = raw_face.vertex_ids;
    auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids)
      vfc = vfc + grid.vertices[fvid];
    newFace.centroid = vfc / double(newFace.vertex_ids.size());

    // A polygon face is just a line so we can just grab
    // the first vertex and form a vector with the face
    // centroid. The normal is then just khat
    // cross-product with this vector.
    uint64_t fvid = newFace.vertex_ids[0];
    auto vec_vvc = grid.vertices[fvid] - newFace.centroid;

    newFace.normal = chi_mesh::Vector3(0.0,0.0,1.0).Cross(vec_vvc);
    newFace.normal.Normalize();

    poly_cell->faces.push_back(newFace);
  }

  grid.cells.push_back(poly_cell);
}

//###################################################################
/**Adds a polyhedron to the grid from a light-weight cell.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  AddPolyhedronToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell &raw_cell,
    const chi_mesh::Cell &temp_cell,
    chi_mesh::MeshContinuum &grid)
{
  auto polyh_cell = new chi_mesh::Cell(CellType::POLYHEDRON, raw_cell.sub_type);
  polyh_cell->centroid = temp_cell.centroid;
  polyh_cell->global_id = temp_cell.global_id;
  polyh_cell->partition_id = temp_cell.partition_id;
  polyh_cell->material_id = temp_cell.material_id;

  polyh_cell->vertex_ids = raw_cell.vertex_ids;

  for (auto& raw_face : raw_cell.faces)
  {
    chi_mesh::CellFace newFace;

    newFace.has_neighbor = raw_face.has_neighbor;
    newFace.neighbor_id = raw_face.neighbor;

    newFace.vertex_ids = raw_face.vertex_ids;
    auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids)
      vfc = vfc + grid.vertices[fvid];
    newFace.centroid = vfc / double(newFace.vertex_ids.size());

    newFace.normal = chi_mesh::Normal(0.0,0.0,0.0);
    size_t last_vert_ind = newFace.vertex_ids.size()-1;
    for (size_t fv=0; fv<newFace.vertex_ids.size(); ++fv)
    {
      uint64_t fvid_m = newFace.vertex_ids[fv];
      uint64_t fvid_p = (fv == last_vert_ind)? newFace.vertex_ids[0] :
                                               newFace.vertex_ids[fv+1];
      auto leg_m = grid.vertices[fvid_m] - newFace.centroid;
      auto leg_p = grid.vertices[fvid_p] - newFace.centroid;

      auto vn = leg_m.Cross(leg_p);

      newFace.normal = newFace.normal + vn.Normalized();
    }
    newFace.normal = (newFace.normal/double(newFace.vertex_ids.size())).Normalized();

    polyh_cell->faces.push_back(newFace);
  }

  grid.cells.push_back(polyh_cell);
}