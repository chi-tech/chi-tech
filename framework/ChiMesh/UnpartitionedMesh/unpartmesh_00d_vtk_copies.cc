#include "chi_unpartitioned_mesh.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>

#include <vtkCellData.h>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Copies the vtk data structures to the current object's internal
 * data.*/
void chi_mesh::UnpartitionedMesh::
  CopyUGridCellsAndPoints(vtkUnstructuredGrid& ugrid,
                          const double scale)
{
  const std::string fname =
    "chi_mesh::UnpartitionedMesh::CopyUGridCellsAndPoints";

  uint64_t total_cell_count  = ugrid.GetNumberOfCells();
  uint64_t total_point_count = ugrid.GetNumberOfPoints();

  chi::log.Log()
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Push cells
  for (size_t c=0; c<total_cell_count; ++c)
  {
    auto vtk_cell = ugrid.GetCell(static_cast<vtkIdType>(c));
    auto vtk_celldim = vtk_cell->GetCellDimension();

    if (vtk_celldim == 3)
      raw_cells_.push_back(CreateCellFromVTKPolyhedron(vtk_cell));
    else if (vtk_celldim == 2)
      raw_cells_.push_back(CreateCellFromVTKPolygon(vtk_cell));
    else if (vtk_celldim == 1)
      raw_cells_.push_back(CreateCellFromVTKLine(vtk_cell));
    else if (vtk_celldim == 0)
      raw_cells_.push_back(CreateCellFromVTKVertex(vtk_cell));
    else
      throw std::logic_error(fname + ": Unsupported cell dimension.");
  }//for c

  //======================================== Push points
  for (size_t p=0; p<total_point_count; ++p)
  {
    auto point = ugrid.GetPoint(static_cast<vtkIdType>(p));

    point[0] = point[0]*scale;
    point[1] = point[1]*scale;
    point[2] = point[2]*scale;

    vertices_.emplace_back(point[0], point[1], point[2]);

    if (not bound_box_)
      bound_box_ = std::shared_ptr<BoundBox>(new BoundBox{point[0],point[0],
                                                          point[1],point[1],
                                                          point[2],point[2]});

    bound_box_->xmin = std::min(bound_box_->xmin, point[0]);
    bound_box_->xmax = std::max(bound_box_->xmax, point[0]);
    bound_box_->ymin = std::min(bound_box_->ymin, point[1]);
    bound_box_->ymax = std::max(bound_box_->ymax, point[1]);
    bound_box_->zmin = std::min(bound_box_->zmin, point[2]);
    bound_box_->zmax = std::max(bound_box_->zmax, point[2]);
  }
}


//###################################################################
/**Set material-ids based on block-wise material ids.*/
void chi_mesh::UnpartitionedMesh::
  SetMaterialIDsFromBlocks(const std::vector<uint64_t> &block_mat_ids)
{
  size_t total_cell_count = raw_cells_.size();
  for (size_t c=0; c<total_cell_count; ++c)
  {
    auto cell = raw_cells_[c];

    int mat_id=-1;
    uint64_t prev_block_lim = 0;
    for (auto block_lim : block_mat_ids)
    {
      ++mat_id;
      if (c<block_lim and c>=prev_block_lim) break;

      prev_block_lim=block_lim;
    }

    cell->material_id = mat_id;
  }
}


//###################################################################
/**Set material-ids from list.*/
void chi_mesh::UnpartitionedMesh::
  SetMaterialIDsFromList(const std::vector<int> &material_ids)
{
  const size_t total_cell_count = raw_cells_.size();

  for (size_t c = 0; c < total_cell_count; ++c)
    raw_cells_[c]->material_id = material_ids[c];
}


//###################################################################
/**Set boundary-ids from boundary grid_blocks.*/
void chi_mesh::UnpartitionedMesh::
  SetBoundaryIDsFromBlocks(
    std::vector<vtkUGridPtrAndName> &bndry_grid_blocks)
{
  const double EPSILON = 1.0e-12;
  //======================================== Build boundary faces
  std::vector<LightWeightFace*> bndry_faces;
  for (auto& cell_ptr : raw_cells_)
    for (auto& face : cell_ptr->faces)
      if (not face.has_neighbor)
        bndry_faces.push_back(&face);

  chi::log.Log() << "Number of boundary faces: " << bndry_faces.size();

  //======================================== Build boundary vertex ids
  std::set<uint64_t> bndry_vids_set;
  for (const auto& face_ptr : bndry_faces)
    for (const auto vid : face_ptr->vertex_ids)
      bndry_vids_set.insert(vid);

  //======================================== Process each boundary block
  uint64_t bndry_id = 0;
  for (const auto& ugrid_name : bndry_grid_blocks)
  {
    auto ugrid = ugrid_name.first;

    mesh_options_.boundary_id_map[bndry_id] = ugrid_name.second;

    //================================= Build vertex map
    bool mapping_failed = false;
    std::vector<size_t> vertex_map(ugrid->GetNumberOfPoints(),0);
    for (size_t p=0; p<ugrid->GetNumberOfPoints(); ++p)
    {
      chi_mesh::Vector3 point;
      ugrid->GetPoint(static_cast<vtkIdType>(p), &point.x);

      bool map_found = false;
      for (const auto vid : bndry_vids_set)
        if ((point - vertices_[vid]).NormSquare() < EPSILON)
        {
          vertex_map[p] = vid;
          map_found = true;
          break;
        }

      if (not map_found)
      {
        chi::log.Log0Warning() << "chi_mesh::UnpartitionedMesh::"
          "SetBoundaryIDsFromBlocks: Failed to map a vertex. " +
          point.PrintStr() + " for boundary " + ugrid_name.second +
          " therefore the boundary assignment was skipped.";
        mapping_failed = true;
        break;
      }
    }//for point in boundary block

    if (mapping_failed) continue;

    //================================= Build vertex subscriptions
    std::map<uint64_t , std::set<size_t>> vertex_face_subs;
    for (size_t f=0; f<bndry_faces.size(); ++f)
      for (const auto vid : bndry_faces[f]->vertex_ids)
        vertex_face_subs[vid].insert(f);

    //================================= Process each cell in bndry block
    size_t num_faces_boundarified = 0;
    auto& bndry_block = ugrid_name.first;
    size_t num_bndry_block_cells = bndry_block->GetNumberOfCells();
    for (size_t bc=0; bc<num_bndry_block_cells; ++bc)
    {
      auto bndry_cell = bndry_block->GetCell(static_cast<vtkIdType>(bc));

      //========================== Build list of face candidates
      //                           and vertex set
      std::set<size_t> face_ids_short_list;
      std::set<uint64_t> bndry_cell_id_set;
      size_t num_points = bndry_cell->GetNumberOfPoints();
      for (size_t p=0; p<num_points; ++p)
      {
        auto point_id = bndry_cell->GetPointId(static_cast<int>(p));
        bndry_cell_id_set.insert(vertex_map[point_id]);
        for (size_t face_id : vertex_face_subs[vertex_map[point_id]])
          face_ids_short_list.insert(face_id);
      }//for point p

      for (size_t face_id : face_ids_short_list)
      {
        auto& face = bndry_faces[face_id];
        const auto& face_vids = face->vertex_ids;
        std::set<uint64_t> face_id_set(face_vids.begin(), face_vids.end());

        if (face_id_set == bndry_cell_id_set)
        {
          face->neighbor = bndry_id;
          ++num_faces_boundarified;
        }
      }//for face_id
    }//for boundary cell bc

    chi::log.Log()
    << "UnpartitionedMesh: assigned " << num_faces_boundarified
    << " to boundary id " << bndry_id
    << " with name " << ugrid_name.second;

    ++bndry_id;
  }//for boundary_block
}
