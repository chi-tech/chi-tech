#include "chi_unpartitioned_mesh.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include <vtkPolyhedron.h>
#include <vtkHexahedron.h>
#include <vtkTetra.h>

#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Creates a raw polyhedron cell from a vtk-polyhedron.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKPolyhedron(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYHEDRON);

  auto vtk_polyh   = vtkPolyhedron::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polyh->GetNumberOfPoints();
  auto num_cfaces  = vtk_polyh->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polyh->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_polyh->GetFace(f);
    auto num_face_points = vtk_face->GetNumberOfPoints();

    face.vertex_ids.reserve(num_face_points);
    auto face_point_ids = vtk_face->GetPointIds();
    for (int p = 0; p < num_face_points; ++p) {
      int point_id = face_point_ids->GetId(p);
      face.vertex_ids.push_back(point_id);
    }

    polyh_cell->faces.push_back(face);
  }

  return polyh_cell;
}

//###################################################################
/**Creates a raw polyhedron cell from a vtk-hexahedron.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKHexahedron(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYHEDRON);

  auto vtk_hex     = vtkHexahedron::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_hex->GetNumberOfPoints();
  auto num_cfaces  = vtk_hex->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_hex->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_hex->GetFace(f);
    auto num_face_points = vtk_face->GetNumberOfPoints();

    face.vertex_ids.reserve(num_face_points);
    auto face_point_ids = vtk_face->GetPointIds();
    for (int p = 0; p < num_face_points; ++p) {
      int point_id = face_point_ids->GetId(p);
      face.vertex_ids.push_back(point_id);
    }

    polyh_cell->faces.push_back(face);
  }

  return polyh_cell;
}

//###################################################################
/**Creates a raw polyhedron cell from a vtk-tetrahedron.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKTetrahedron(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYHEDRON);

  auto vtk_tet     = vtkTetra::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_tet->GetNumberOfPoints();
  auto num_cfaces  = vtk_tet->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_tet->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_tet->GetFace(f);
    auto num_face_points = vtk_face->GetNumberOfPoints();

    face.vertex_ids.reserve(num_face_points);
    auto face_point_ids = vtk_face->GetPointIds();
    for (int p = 0; p < num_face_points; ++p) {
      int point_id = face_point_ids->GetId(p);
      face.vertex_ids.push_back(point_id);
    }

    polyh_cell->faces.push_back(face);
  }

  return polyh_cell;
}

//###################################################################
/**Creates a raw polyhedron cell from a vtk-polygon.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKPolygon(vtkCell *vtk_cell)
{
  auto poly_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_polygon    = vtkPolygon::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polygon->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polygon->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    poly_cell->vertex_ids.push_back(point_id);
  }//for p

  poly_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v0_id = poly_cell->vertex_ids[f];
    auto v1_id = (f<(num_cfaces-1))? poly_cell->vertex_ids[f+1] :
                                     poly_cell->vertex_ids[0];

    face.vertex_ids.reserve(2);
    face.vertex_ids.push_back(v0_id);
    face.vertex_ids.push_back(v1_id);

    poly_cell->faces.push_back(face);
  }

  return poly_cell;
}

//###################################################################
/**Creates a raw polyhedron cell from a vtk-quad.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKQuad(vtkCell *vtk_cell)
{
  auto poly_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_polygon    = vtkPolygon::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polygon->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polygon->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    poly_cell->vertex_ids.push_back(point_id);
  }//for p

  poly_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v0_id = poly_cell->vertex_ids[f];
    auto v1_id = (f<(num_cfaces-1))? poly_cell->vertex_ids[f+1] :
                 poly_cell->vertex_ids[0];

    face.vertex_ids.reserve(2);
    face.vertex_ids.push_back(v0_id);
    face.vertex_ids.push_back(v1_id);

    poly_cell->faces.push_back(face);
  }

  return poly_cell;
}

//###################################################################
/**Creates a raw polyhedron cell from a vtk-triangle.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKTriangle(vtkCell *vtk_cell)
{
  auto poly_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_polygon    = vtkPolygon::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polygon->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polygon->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    poly_cell->vertex_ids.push_back(point_id);
  }//for p

  poly_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v0_id = poly_cell->vertex_ids[f];
    auto v1_id = (f<(num_cfaces-1))? poly_cell->vertex_ids[f+1] :
                 poly_cell->vertex_ids[0];

    face.vertex_ids.reserve(2);
    face.vertex_ids.push_back(v0_id);
    face.vertex_ids.push_back(v1_id);

    poly_cell->faces.push_back(face);
  }

  return poly_cell;
}

//###################################################################
/**Establishes neighbor connectivity for the light-weight mesh.*/
void chi_mesh::UnpartitionedMesh::BuildMeshConnectivity()
{
  //======================================== Populate vertex
  //                                                   subscriptionns
  std::vector<std::set<size_t>> vertex_subs(vertices.size());
  uint64_t cur_cell_id=0;
  for (auto& cell : raw_cells)
  {
    for (auto vid : cell->vertex_ids)
      vertex_subs[vid].insert(cur_cell_id);
    ++cur_cell_id;
  }

  int num_bndry_faces = 0;
  int cell_cnt = 0;
  for (auto& cell : raw_cells)
  {
    for (auto& face : cell->faces)
    {
      if (face.neighbor < 0) ++num_bndry_faces;
      face.neighbor = -1;
    }
  }

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "before connectivity: " << num_bndry_faces;

  //======================================== Establish connectivity
  chi_log.Log() << "Establishing cell connectivity.";
  std::set<size_t> cells_to_search;
  cur_cell_id=0;
  for (auto& cell : raw_cells)
  {
    cells_to_search.clear();
    for (uint64_t vid : cell->vertex_ids)
      for (uint64_t cell_id : vertex_subs[vid])
        if (cell_id != cur_cell_id)
          cells_to_search.insert(cell_id);

    for (auto& cur_cell_face : cell->faces)
    {
      if (cur_cell_face.neighbor >= 0 ) continue;

      std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                cur_cell_face.vertex_ids.end());

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto adj_cell = raw_cells[adj_cell_id];

        for (auto& adj_cell_face : adj_cell->faces)
        {
          if (adj_cell_face.neighbor >= 0) continue;
          std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                    adj_cell_face.vertex_ids.end());

          if (cfvids == afvids)
          {
            cur_cell_face.neighbor = adj_cell_id;
            adj_cell_face.neighbor = cur_cell_id;
            goto face_neighbor_found;
          }
        }//for adjacent cell face
      }
      face_neighbor_found:;
    }//for face

    ++cur_cell_id;
  }//for cell

  chi_log.Log() << "Done establishing cell connectivity.";

  num_bndry_faces = 0;
  for (auto cell : raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "after connectivity: " << num_bndry_faces;
}

/**Compute centroids for all cells.*/
void chi_mesh::UnpartitionedMesh::ComputeCentroidsAndCheckQuality()
{
  chi_log.Log() << "Computing cell-centroids.";
  for (auto cell : raw_cells)
  {
    cell->centroid = chi_mesh::Vertex(0.0,0.0,0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid += *vertices[vid];

    cell->centroid = cell->centroid/(cell->vertex_ids.size());
  }
  chi_log.Log() << "Done computing cell-centroids.";

  chi_log.Log() << "Checking cell-center-to-face orientations";
  size_t check0_corrections=0;
  for (auto cell : raw_cells)
    if (cell->type == CellType::POLYHEDRON)
      for (auto& face : cell->faces)
      {
        chi_mesh::Vector3 face_centroid;
        for (uint64_t vid : face.vertex_ids)
          face_centroid += *vertices[vid];
        face_centroid /= face.vertex_ids.size();

        if (face.vertex_ids.size()<2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face with less than"
                                 " 2 vertices on a face, making normal computation impossible.");

        const auto& fv1 = *vertices[face.vertex_ids[0]];
        const auto& fv2 = *vertices[face.vertex_ids[1]];

        auto E0 = fv1-face_centroid;
        auto E1 = fv2-face_centroid;
        auto n  = E0.Cross(E1); n.Normalize();

        auto CC = face_centroid - cell->centroid;

        if (n.Dot(CC) < 0.0)
        {
          std::vector<uint64_t> reversed_verts(face.vertex_ids.rbegin(),
                                               face.vertex_ids.rend());
          face.vertex_ids = reversed_verts;
          ++check0_corrections;
        }
      }//for face

  if (check0_corrections > 0)
    chi_log.Log(LOG_ALLWARNING)
      << "Cell-center-to-face orientation detected " << check0_corrections
      << " faces that violated quality requirements. An attempt to fix it"
      << " was made.";
  chi_log.Log() << "Done checking cell-center-to-face orientations";
}