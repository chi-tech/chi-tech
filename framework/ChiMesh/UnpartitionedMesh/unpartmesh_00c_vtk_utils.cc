#include "chi_unpartitioned_mesh.h"

#include <vtkPolyhedron.h>
#include <vtkHexahedron.h>
#include <vtkTetra.h>

#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>

#include <vtkLine.h>

#include <vtkVertex.h>

//###################################################################
/**Creates a raw polyhedron cell from a vtk-polyhedron.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKPolyhedron(vtkCell *vtk_cell)
{
  auto polyh_cell  = new LightWeightCell(CellType::POLYHEDRON,
                                         CellType::POLYHEDRON);

  auto vtk_polyh   = vtkPolyhedron::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polyh->GetNumberOfPoints();
  auto num_cfaces  = vtk_polyh->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polyh->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
      uint64_t point_id = face_point_ids->GetId(p);
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
  auto polyh_cell  = new LightWeightCell(CellType::POLYHEDRON,
                                         CellType::HEXAHEDRON);

  auto vtk_hex     = vtkHexahedron::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_hex->GetNumberOfPoints();
  auto num_cfaces  = vtk_hex->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_hex->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
      uint64_t point_id = face_point_ids->GetId(p);
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
  auto polyh_cell  = new LightWeightCell(CellType::POLYHEDRON,
                                         CellType::TETRAHEDRON);

  auto vtk_tet     = vtkTetra::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_tet->GetNumberOfPoints();
  auto num_cfaces  = vtk_tet->GetNumberOfFaces();

  polyh_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_tet->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
      uint64_t point_id = face_point_ids->GetId(p);
      face.vertex_ids.push_back(point_id);
    }

    polyh_cell->faces.push_back(face);
  }

  return polyh_cell;
}

//###################################################################
/**Creates a raw polygon cell from a vtk-polygon.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKPolygon(vtkCell *vtk_cell)
{
  auto poly_cell   = new LightWeightCell(CellType::POLYGON,
                                         CellType::POLYGON);

  auto vtk_polygon = vtkPolygon::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_polygon->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_polygon->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
/**Creates a raw polygon cell from a vtk-quad.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKQuad(vtkCell *vtk_cell)
{
  auto poly_cell   = new LightWeightCell(CellType::POLYGON,
                                         CellType::QUADRILATERAL);

  auto vtk_quad    = vtkQuad::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_quad->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_quad->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
/**Creates a raw polygon cell from a vtk-triangle.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKTriangle(vtkCell *vtk_cell)
{
  auto poly_cell   = new LightWeightCell(CellType::POLYGON,
                                         CellType::TRIANGLE);

  auto vtk_triangle= vtkTriangle::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_triangle->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  poly_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_triangle->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
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
/**Creates a raw slab cell from a vtk-line.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKLine(vtkCell *vtk_cell)
{
  auto slab_cell   = new LightWeightCell(CellType::SLAB,
                                         CellType::SLAB);

  auto vtk_line    = vtkLine::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_line->GetNumberOfPoints();
  auto num_cfaces  = num_cpoints;

  slab_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_line->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    slab_cell->vertex_ids.push_back(point_id);
  }//for p

  slab_cell->faces.reserve(num_cfaces);
  for (int f=0; f<num_cfaces; ++f)
  {
    LightWeightFace face;

    auto v_id = slab_cell->vertex_ids[f];

    face.vertex_ids.reserve(1);
    face.vertex_ids.push_back(v_id);

    slab_cell->faces.push_back(face);
  }

  return slab_cell;
}

//###################################################################
/**Creates a raw point cell from a vtk-vertex.*/
chi_mesh::UnpartitionedMesh::LightWeightCell* chi_mesh::UnpartitionedMesh::
  CreateCellFromVTKVertex(vtkCell *vtk_cell)
{
  auto point_cell  = new LightWeightCell(CellType::GHOST,
                                         CellType::POINT);

  auto vtk_vertex  = vtkVertex::SafeDownCast(vtk_cell);
  auto num_cpoints = vtk_vertex->GetNumberOfPoints();

  point_cell->vertex_ids.reserve(num_cpoints);
  auto point_ids   = vtk_vertex->GetPointIds();
  for (int p=0; p<num_cpoints; ++p)
  {
    uint64_t point_id = point_ids->GetId(p);
    point_cell->vertex_ids.push_back(point_id);
  }//for p

  return point_cell;
}