#ifndef CHITECH_LBS_SOLVER_H
#define CHITECH_LBS_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/LinearSolver/linear_solver.h"
#include "lbs_structs.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"

#include "A_LBSSolver/PointSource/lbs_point_source.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_group.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

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

  LBSSolver (const LBSSolver&) = delete;
  LBSSolver& operator= (const LBSSolver&) = delete;

  double LastRestartWrite() const;
  double& LastRestartWrite();

  lbs::Options& Options();

  size_t NumMoments() const;
  size_t NumGroups() const;

  void AddGroup(int id);
  const std::vector<LBSGroup>& Groups() const;

  void AddGroupset();
  std::vector<LBSGroupset>& Groupsets();
  const std::vector<LBSGroupset>& Groupsets() const;

  const chi_math::SpatialDiscretization& SpatialDiscretization() const;

  void AddPointSource(PointSource psrc);
  void ClearPointSources();
  const std::vector<PointSource>& PointSources() const;

  std::map<uint64_t, BoundaryPreference>& BoundaryPreferences();

  const chi_math::UnknownManager& UnknownManager() const;

  size_t LocalNodeCount() const;
  size_t GlobalNodeCount() const;

  std::vector<double>& QMomentsLocal();
  std::vector<double>& ExtSrcMomentsLocal();
  std::vector<double>& PhiOldLocal();
  std::vector<double>& PhiNewLocal();

  //01
  void Initialize() override;
protected:
  //01a
  virtual void PerformInputChecks();
  //01b
  void PrintSimHeader();
public:
  //01c
  void InitMaterials();
protected:
  //01d
  virtual void InitializeSpatialDiscretization();
  void ComputeUnitIntegrals();
  //01e
  void InitializeGroupsets();
  //01f
  void ComputeNumberOfMoments();
  //01g
  virtual void InitializeParrays();
  //01h
  void InitializeBoundaries();
public:
  //01i
  void InitializePointSources();

//03d
  void InitWGDSA(LBSGroupset& groupset);
public:
  void ExecuteWGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& ref_phi_old,
                                   const std::vector<double>& ref_phi_new,
                                   std::vector<double>& delta_phi_local);

  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);

  void DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);
protected:
  static void CleanUpWGDSA(LBSGroupset& groupset);

  //03e
  void InitTGDSA(LBSGroupset& groupset);
public:
  void ExecuteTGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& ref_phi_old,
                                   const std::vector<double>& ref_phi_new,
                                   std::vector<double>& delta_phi_local);
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);
  void DisAssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);
protected:
  static void CleanUpTGDSA(LBSGroupset& groupset);
  //04 File IO
public:
  //04a
  void WriteRestartData(const std::string& folder_name,
                        const std::string& file_base);
  void ReadRestartData(const std::string& folder_name, const std::string& file_base);
  //04b
  void WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                                  const std::string& file_base);
  void ReadGroupsetAngularFluxes(LBSGroupset& groupset,
                                 const std::string& file_base);

  //04c
  std::vector<double> MakeSourceMomentsFromPhi();
  void WriteFluxMoments(const std::string& file_base,
                        const std::vector<double>& flux_moments);
  void ReadFluxMoments(const std::string& file_base,
                       std::vector<double>& flux_moments,
                       bool single_file=false);

  //05a
  void UpdateFieldFunctions();

  //Vector assembly
  virtual void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                                 PhiSTLOption which_phi);

  virtual void SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                                 PhiSTLOption which_phi);

  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             const std::vector<double>& x_src,
                                             std::vector<double>& y);

  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             PhiSTLOption from_which_phi,
                                             PhiSTLOption to_which_phi);

  virtual void SetGroupScopedPETScVecFromPrimarySTLvector(
    int first_group_id,
    int last_group_id, Vec x,
    const std::vector<double>& y);

  virtual void SetPrimarySTLvectorFromGroupScopedPETScVec(
    int first_group_id,
    int last_group_id, Vec x_src,
    std::vector<double>& y);
};


}//namespace lbs

#endif //CHITECH_LBS_SOLVER_H
