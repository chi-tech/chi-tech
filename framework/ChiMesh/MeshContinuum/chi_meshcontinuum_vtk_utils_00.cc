#include "chi_meshcontinuum_vtk_utils.h"

#include "chi_meshcontinuum.h"

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

//###################################################################
/**Uploads vertices and cells to an unstructured grid.*/
void chi_mesh::UploadCellGeometry(const chi_mesh::MeshContinuum &grid,
                                  const chi_mesh::Cell &cell,
                                  int64_t &node_counter,
                                  vtkNew<vtkPoints> &points,
                                  vtkNew<vtkUnstructuredGrid> &ugrid)
{
  size_t num_verts = cell.vertex_ids.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v=0; v<num_verts; v++)
  {
    uint64_t vgi = cell.vertex_ids[v];
    std::vector<double> d_node(3);
    d_node[0] = grid.vertices[vgi].x;
    d_node[1] = grid.vertices[vgi].y;
    d_node[2] = grid.vertices[vgi].z;

    points->InsertPoint(node_counter,d_node.data());
    cell_vids[v] = node_counter++;
  }

  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    ugrid->InsertNextCell(VTK_LINE,
                          static_cast<vtkIdType>(num_verts),
                          cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    ugrid->InsertNextCell(VTK_POLYGON,
                          static_cast<vtkIdType>(num_verts),
                          cell_vids.data());
  }
  if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    // Build polyhedron faces
    std::vector<vtkIdType> faces_vids;

    size_t num_faces = cell.faces.size();
    for (auto& face : cell.faces)
    {
      size_t num_fverts = face.vertex_ids.size();
      std::vector<vtkIdType> face_info(num_fverts);
      for (size_t fv=0; fv<num_fverts; fv++)
      {
        size_t v = 0;
        for (size_t cv=0; cv<num_verts; ++cv)
          if (cell.vertex_ids[cv] == face.vertex_ids[fv])
          { v = cv; break; }

        face_info[fv] = cell_vids[v];
      }

      faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
      for (auto vid : face_info)
        faces_vids.push_back(vid);
    }//for f

    ugrid->InsertNextCell(VTK_POLYHEDRON,
                          static_cast<vtkIdType>(num_verts),
                          cell_vids.data(),
                          static_cast<vtkIdType>(num_faces),
                          faces_vids.data());
  }//polyhedron
}