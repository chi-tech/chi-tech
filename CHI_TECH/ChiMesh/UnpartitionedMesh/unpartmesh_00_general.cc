#include "chi_unpartitioned_mesh.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include <vtkPolyhedron.h>
#include <vtkHexahedron.h>
#include <vtkTetra.h>

#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>

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