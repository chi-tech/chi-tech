#ifndef CHITECH_CHI_MESHCONTINUUM_VTK_UTILS_H
#define CHITECH_CHI_MESHCONTINUUM_VTK_UTILS_H

#include <cstdint>
#include <string>

template<class T>
class vtkNew;
class vtkPoints;
class vtkUnstructuredGrid;

namespace chi_mesh
{
  class MeshContinuum;
  class Cell;

  void UploadCellGeometry(const chi_mesh::MeshContinuum& grid,
                          const chi_mesh::Cell &cell,
                          int64_t& node_counter,
                          vtkNew<vtkPoints>& points,
                          vtkNew<vtkUnstructuredGrid>& ugrid);

  vtkNew<vtkUnstructuredGrid>
    PrepareVtkUnstructuredGrid(const chi_mesh::MeshContinuum& grid);

  void WritePVTUFiles(vtkNew<vtkUnstructuredGrid> &ugrid,
                      const std::string &file_base_name);

}//namespace chi_mesh

#endif //CHITECH_CHI_MESHCONTINUUM_VTK_UTILS_H
