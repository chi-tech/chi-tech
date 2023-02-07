#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "lbs_structs.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"

#include "LBSSteadyState/PointSource/lbs_point_source.h"
#include "LBSSteadyState/Groupset/lbs_group.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

#include <petscksp.h>

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;
//typedef sweep_namespace::SweepScheduler MainSweepScheduler;

namespace lbs
{
//################################################################### Class def
/**A neutral particle transport solver.*/
class SteadyStateSolver : public chi_physics::Solver
{
protected:
  typedef chi_mesh::sweep_management::CellFaceNodalMapping CellFaceNodalMapping;
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

 public:
  //00
  explicit SteadyStateSolver(const std::string& in_text_name);
  ~SteadyStateSolver() override =default;

  SteadyStateSolver (const SteadyStateSolver&) = delete;
  SteadyStateSolver& operator= (const SteadyStateSolver&) = delete;

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





public:
  //02
  void Execute() override;
protected:
  virtual void SolveGroupset(LBSGroupset& groupset);

  //03a
  void ComputeSweepOrderings(LBSGroupset& groupset) const;
  //03aa
  static
  std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(const chi_mesh::MeshContinuum& grid,
                                 const chi_math::AngularQuadrature& quadrature,
                                 AngleAggregationType agg_type,
                                 lbs::GeometryType lbs_geo_type);
  //03b
  void InitFluxDataStructures(LBSGroupset& groupset);

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
  //03f
  void ResetSweepOrderings(LBSGroupset& groupset);

  //04 File IO
  //04a
public:
  void WriteRestartData(std::string folder_name, std::string file_base);
  void ReadRestartData(std::string folder_name, std::string file_base);

public:
  //05
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

  //Iterative Operations
  void SetSource(LBSGroupset& groupset,
                 std::vector<double>& destination_q,
                 const std::vector<double>& phi,
                 SourceFlags source_flags);
protected:
  double ComputePiecewiseChange(LBSGroupset& groupset);
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);
  double ComputeFissionProduction(const std::vector<double>& phi);
public:
  virtual double ComputeFissionRate(bool previous);
protected:
  //Iterative Methods
  bool ClassicRichardson(LBSGroupset& groupset,
                         chi_mesh::sweep_management::SweepScheduler& sweep_scheduler,
                         SourceFlags source_flags,
                         const SetSourceFunction& set_source_function,
                         bool log_info = true);
  bool Krylov(LBSGroupset& groupset,
              chi_mesh::sweep_management::SweepScheduler& sweep_scheduler,
              SourceFlags lhs_src_scope,
              SourceFlags rhs_src_scope,
              const SetSourceFunction& set_source_function,
              bool log_info = true);

  //Vector assembly
public:
  void SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x,
                                         const std::vector<double>& y,
                                         bool with_delayed_psi=false);
  void SetGSSTLvectorFromPrimarySTLvector(LBSGroupset& groupset,
                                          std::vector<double>& x,
                                          const std::vector<double>& y,
                                          bool with_delayed_psi=false);
  void SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src,
                                         std::vector<double>& y,
                                         bool with_delayed_psi=false);
  void SetPrimarySTLvectorFromGSSTLvector(LBSGroupset& groupset,
                                          const std::vector<double>& x_src,
                                          std::vector<double>& y,
                                          bool with_delayed_psi=false);
  void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                     const std::vector<double>& x_src,
                                     std::vector<double>& y,
                                     bool with_delayed_psi=false);
protected:
  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
public:
  void ComputeBalance();

protected:
  //precursors
  void ComputePrecursors();
};

}

#endif
