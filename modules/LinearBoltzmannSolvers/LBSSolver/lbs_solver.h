#ifndef CHITECH_LBS_SOLVER_H
#define CHITECH_LBS_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/LinearSolver/linear_solver.h"
#include "lbs_structs.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"

#include "LBSSteadyState/PointSource/lbs_point_source.h"
#include "LBSSteadyState/Groupset/lbs_group.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

#include <petscksp.h>

namespace lbs
{
  template<class MatType, class VecType, class SolverType>
  class AGSLinearSolver;
  template<class MatType, class VecType, class SolverType>
  class WGSLinearSolver;
}

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;

namespace lbs
{

//################################################################### Class def
/**Base class for all Linear Boltzmann Solvers.*/
class LBSSolver : public chi_physics::Solver
{
protected:
  typedef chi_mesh::sweep_management::CellFaceNodalMapping CellFaceNodalMapping;
  typedef std::shared_ptr<AGSLinearSolver<Mat,Vec,KSP>> AGSLinSolverPtr;
  typedef std::shared_ptr<chi_math::LinearSolver<Mat,Vec,KSP>> LinSolvePtr;

  size_t source_event_tag_=0;
  double last_restart_write_=0.0;

  lbs::Options options_;
  size_t num_moments_ = 0;
  size_t num_groups_ = 0;
  size_t num_precursors_ = 0;
  size_t max_precursors_per_material_ = 0;

  std::vector<LBSGroup> groups_;
  std::vector<LBSGroupset> groupsets_;
  std::vector<PointSource> point_sources_;

  std::map<int,XSPtr>           matid_to_xs_map_;
  std::map<int,IsotropicSrcPtr> matid_to_src_map_;

  std::shared_ptr<chi_math::SpatialDiscretization> discretization_ = nullptr;
  chi_mesh::MeshContinuumPtr grid_ptr_;
  std::vector<CellFaceNodalMapping> grid_nodal_mappings_;
  std::vector<UnitCellMatrices> unit_cell_matrices_;
  std::vector<lbs::CellLBSView> cell_transport_views_;

  std::map<uint64_t, BoundaryPreference>           boundary_preferences_;
  std::map<uint64_t, std::shared_ptr<SweepBndry>>  sweep_boundaries_;

  chi_math::UnknownManager flux_moments_uk_man_;

  size_t max_cell_dof_count_ = 0;
  uint64_t local_node_count_ = 0;
  uint64_t glob_node_count_ = 0;

  std::vector<double> q_moments_local_, ext_src_moments_local_;
  std::vector<double> phi_new_local_, phi_old_local_;
  std::vector<std::vector<double>> psi_new_local_;
  std::vector<double> precursor_new_local_;

  SetSourceFunction active_set_source_function_;

  std::vector<AGSLinSolverPtr> ags_solvers_;
  std::vector<LinSolvePtr>     wgs_solvers_;
  AGSLinSolverPtr              primary_ags_solver_;
public:
  explicit LBSSolver(const std::string& text_name);
};


}//namespace lbs

#endif //CHITECH_LBS_SOLVER_H
