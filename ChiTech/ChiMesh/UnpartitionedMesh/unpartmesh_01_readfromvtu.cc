#include "chi_unpartitioned_mesh.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <fstream>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

#include <vtkPolyhedron.h>
#include <vtkHexahedron.h>
#include <vtkTetra.h>

#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>

#include <vtkLine.h>

#include <vtkVertex.h>

#include <vtkCleanUnstructuredGrid.h>

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

//###################################################################
/**Reads a VTK unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromVTU(const chi_mesh::UnpartitionedMesh::Options &options)
{
  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);

  if (!file.is_open())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to open file: "<< options.file_name <<" in call "
      << "to ReadFromVTU \n";
    exit(EXIT_FAILURE);
  }
  file.close();

  //======================================== Read the file
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());
  chi_log.Log(LOG_0)
    << "Reading VTU file     : \""
    << options.file_name << "\".";

  reader->Update();

  chi_log.Log(LOG_0)
    << "Done reading VTU file: \""
    << options.file_name << "\".";

  //======================================== Remove duplicate vertices
  auto cleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
  cleaner->SetInputData(reader->GetOutput());
  cleaner->Update();
  auto ugrid = cleaner->GetOutput();

  auto total_cell_count  = ugrid->GetNumberOfCells();
  auto total_point_count = ugrid->GetNumberOfPoints();

  chi_log.Log(LOG_0)
    << "Clean grid num cells and points: "
    << total_cell_count << " "
    << total_point_count;

  //======================================== Determine if reading
  //                                         cell identifiers
  if (options.material_id_fieldname != options.boundary_id_fieldname)
  {
    chi_log.Log(LOG_ALLERROR)
      << "The VTU reader expects material identifiers and boundary identifiers "
      << "to be defined in the same field.";
    std::exit(EXIT_FAILURE);
  }

  vtkDataArray* cell_id_array_ptr = nullptr;
  if (options.material_id_fieldname.empty())
  {
    chi_log.Log(LOG_0)
      << "A user-supplied field name from which to recover cell identifiers "
      << "has not been provided. Only the mesh will be read.";
  }
  else
  {
    chi_log.Log(LOG_0)
      << "A user-supplied field name from which to recover cell identifiers "
      << "has been provided. The mesh will be read and both material ID and "
      << "boundary ID will be read from the vtkCellData field with name : \""
      << options.material_id_fieldname << "\".";

    const auto vtk_abstract_array_ptr =
      ugrid->GetCellData()->GetAbstractArray(options.material_id_fieldname.c_str());
    if (!vtk_abstract_array_ptr)
    {
      chi_log.Log(LOG_ALLERROR)
        << "The VTU file : \"" << options.file_name << "\" "
        << "does not contain a vtkCellData field of name : \""
        << options.material_id_fieldname << "\".";
      std::exit(EXIT_FAILURE);
    }

    cell_id_array_ptr = vtkArrayDownCast<vtkDataArray>(vtk_abstract_array_ptr);
    if (!cell_id_array_ptr)
    {
      chi_log.Log(LOG_ALLERROR)
        << "The VTU file : \"" << options.file_name << "\" "
        << "with vtkCellData field of name : \""
        << options.material_id_fieldname << "\" "
        << "cannot be downcast to vtkDataArray";
      std::exit(EXIT_FAILURE);
    }

    const auto cell_id_n_tup = cell_id_array_ptr->GetNumberOfTuples();
    if (cell_id_n_tup != total_cell_count)
    {
      chi_log.Log(LOG_ALLERROR)
        << "The VTU file : \"" << options.file_name << "\" "
        << "with vtkCellData field of name : \""
        << options.material_id_fieldname << "\" has n. tuples : "
        << cell_id_n_tup << ", but differs from the value expected : "
        << total_cell_count << ".";
      std::exit(EXIT_FAILURE);
    }

    const auto cell_id_n_val = cell_id_array_ptr->GetNumberOfValues();
    if (cell_id_n_val != total_cell_count)
    {
      chi_log.Log(LOG_ALLERROR)
        << "The VTU file : \"" << options.file_name << "\" "
        << "with vtkCellData field of name : \""
        << options.material_id_fieldname << "\" has n. values : "
        << cell_id_n_val << ", but differs from the value expected : "
        << total_cell_count << ".";
      std::exit(EXIT_FAILURE);
    }
  }

  //======================================== Scan for mesh dimension
  int mesh_dim = 0;
  for (unsigned int c = 0; c < total_cell_count; ++c)
  {
    const auto vtk_cell = ugrid->GetCell(c);
    const auto vtk_cell_dim = vtk_cell->GetCellDimension();
    mesh_dim = std::max(vtk_cell_dim, mesh_dim);
  }

  if (mesh_dim < 1 || mesh_dim > 3)
  {
    chi_log.Log(LOG_ALLERROR)
      << "The VTU file : \"" << options.file_name << "\" "
      << "does not identify a mesh of valid dimension.";
    std::exit(EXIT_FAILURE);
  }

  //======================================== Push cells
  size_t num_polyhedrons  = 0;
  size_t num_hexahedrons  = 0;
  size_t num_tetrahedrons = 0;
  size_t num_polygons     = 0;
  size_t num_quads        = 0;
  size_t num_triangles    = 0;
  size_t num_lines        = 0;
  size_t num_vertices     = 0;
  for (unsigned int c = 0; c < total_cell_count; ++c)
  {
    const auto vtk_cell = ugrid->GetCell(c);
    const auto vtk_celltype = vtk_cell->GetCellType();
    const auto vtk_cell_dim = vtk_cell->GetCellDimension();

    //  construct cell
    LightWeightCell* chi_lwc = nullptr;
    if (vtk_celltype == VTK_POLYHEDRON)
    {
      chi_lwc = CreateCellFromVTKPolyhedron(vtk_cell);
      ++num_polyhedrons;
    }
    else if (vtk_celltype == VTK_HEXAHEDRON)
    {
      chi_lwc = CreateCellFromVTKHexahedron(vtk_cell);
      ++num_hexahedrons;
    }
    else if (vtk_celltype == VTK_TETRA)
    {
      chi_lwc = CreateCellFromVTKTetrahedron(vtk_cell);
      ++num_tetrahedrons;
    }
    else if (vtk_celltype == VTK_POLYGON)
    {
      chi_lwc = CreateCellFromVTKPolygon(vtk_cell);
      ++num_polygons;
    }
    else if (vtk_celltype == VTK_QUAD)
    {
      chi_lwc = CreateCellFromVTKQuad(vtk_cell);
      ++num_quads;
    }
    else if (vtk_celltype == VTK_TRIANGLE)
    {
      chi_lwc = CreateCellFromVTKTriangle(vtk_cell);
      ++num_triangles;
    }
    else if (vtk_celltype == VTK_LINE)
    {
      chi_lwc = CreateCellFromVTKLine(vtk_cell);
      ++num_lines;
    }
    else if (vtk_celltype == VTK_VERTEX)
    {
      chi_lwc = CreateCellFromVTKVertex(vtk_cell);
      ++num_vertices;
    }
    else
      throw std::invalid_argument(std::string(__FUNCTION__) +
                                  ": Unsupported cell type.");

    //  append to appropriate collection
    if (vtk_cell_dim == mesh_dim)
      raw_cells.emplace_back(chi_lwc);
    else if (vtk_cell_dim == mesh_dim-1)
      raw_boundary_cells.emplace_back(chi_lwc);

    //  apply cell identifier
    if (cell_id_array_ptr)
    {
      std::vector<double> cell_id_vec(1);
      cell_id_array_ptr->GetTuple(c, cell_id_vec.data());
      const auto cell_id = (int)cell_id_vec.front();

      if (vtk_cell_dim == mesh_dim)
        raw_cells.back()->material_id = cell_id;
      else if (vtk_cell_dim == mesh_dim-1)
        raw_boundary_cells.back()->material_id = cell_id;
    }
  }//for c
  chi_log.Log() << "Number cells read: " << total_cell_count << "\n"
    << "polyhedrons  : " << num_polyhedrons  << "\n"
    << "hexahedrons  : " << num_hexahedrons  << "\n"
    << "tetrahedrons : " << num_tetrahedrons << "\n"
    << "polygons     : " << num_polygons     << "\n"
    << "quads        : " << num_quads        << "\n"
    << "triangles    : " << num_triangles    << "\n"
    << "lines        : " << num_lines        << "\n"
    << "vertices     : " << num_vertices;

  //======================================== Push points
  for (int p=0; p<total_point_count; ++p)
  {
    auto point = ugrid->GetPoint(p);
    vertices.emplace_back(point[0],point[1],point[2]);

    if (point[0] < bound_box.xmin) bound_box.xmin = point[0];
    if (point[0] > bound_box.xmax) bound_box.xmax = point[0];
    if (point[1] < bound_box.ymin) bound_box.ymin = point[1];
    if (point[1] > bound_box.ymax) bound_box.ymax = point[1];
    if (point[2] < bound_box.zmin) bound_box.zmin = point[2];
    if (point[2] > bound_box.zmax) bound_box.zmax = point[2];
  }

  //======================================== Always do this
  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();
}

