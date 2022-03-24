#include "lbsadj_solver.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

//###################################################################
/**Subscribes cells to QOIs.*/
size_t lbs_adjoint::AdjointSolver::
  SetQOI(const std::string& qoi_name,
         const chi_mesh::LogicalVolume& logical_volume,
         const std::string& lua_function_name)
{
  // Make the designation
  QOIDesignation qoi_designation(qoi_name, logical_volume, lua_function_name);
  // Make empty subscriber list (will be populated during initialize)
  std::vector<size_t> qoi_cell_subscriptions;

  QOI_cell_subscriptions.emplace_back(qoi_designation, qoi_cell_subscriptions);

  return QOI_cell_subscriptions.size() - 1;
}