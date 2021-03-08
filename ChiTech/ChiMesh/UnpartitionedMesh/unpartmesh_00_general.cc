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
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_polygon    = vtkPolygon::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polygon->GetNumberOfPoints();
  auto num_cfaces  = vtk_polygon->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polygon->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_polygon->GetFace(f);
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
/**Creates a raw polyhedron cell from a vtk-quad.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKQuad(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_quad    = vtkQuad::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_quad->GetNumberOfPoints();
  auto num_cfaces  = vtk_quad->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_quad->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_quad->GetFace(f);
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
/**Creates a raw polyhedron cell from a vtk-triangle.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKTriangle(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(chi_mesh::CellType::POLYGON);

  auto vtk_triangle = vtkTriangle::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_triangle->GetNumberOfPoints();
  auto num_cfaces  = vtk_triangle->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_triangle->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    int point_id = point_ids->GetId(p);
    polyh_cell->vertex_ids.push_back(point_id);
  }//for p

  polyh_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;
    auto vtk_face = vtk_triangle->GetFace(f);
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
  for (auto cell : raw_cells)
    for (auto& face : cell->faces)
    {
      if (face.neighbor < 0) ++num_bndry_faces;
      face.neighbor = -1;
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
void chi_mesh::UnpartitionedMesh::ComputeCentroids()
{
  for (auto cell : raw_cells)
  {
    cell->centroid = chi_mesh::Vertex(0.0,0.0,0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid = cell->centroid + *vertices[vid];

    cell->centroid = cell->centroid/(cell->vertex_ids.size());
  }
}