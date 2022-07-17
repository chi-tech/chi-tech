#include "lbsadj_solver.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

//###################################################################
/**Subscribes cells to QOIs.*/
size_t lbs_adjoint::AdjointSolver::
  AddResponseFunction(const std::string& qoi_name,
                      std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
                      const std::string& lua_function_name)
{
  // Make the designation
  ResponseFunctionDesignation qoi_designation(qoi_name, logical_volume, lua_function_name);
  // Make empty subscriber list (will be populated during initialize)
  std::vector<size_t> cell_rf_subscriptions;

  response_functions.emplace_back(qoi_designation, cell_rf_subscriptions);

  return response_functions.size() - 1;
}