#include "chi_meshcontinuum.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_grid_vtk_utils.h"

#include <vtkUnstructuredGrid.h>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Exports just the mesh to VTK format.*/
void chi_mesh::MeshContinuum::ExportCellsToVTK(const std::string& file_base_name) const
{
  chi::log.Log() << "Exporting mesh to VTK files with base " << file_base_name;

  const auto& grid = *this;

  auto ugrid = chi_mesh::PrepareVtkUnstructuredGrid(grid, false);

  chi_mesh::WritePVTUFiles(ugrid, file_base_name);

  chi::log.Log() << "Done exporting mesh to VTK.";
}