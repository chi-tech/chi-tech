#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "Groupset/lbs_groupset.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "lbs_structs.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "PointSource/lbs_point_source.h"
#include "Groupset/lbs_group.h"
#include "Groupset/lbs_groupset.h"

#include <petscksp.h>

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

namespace lbs
{
enum class BoundaryType
{
  VACUUM = 1,
  INCIDENT_ISOTROPIC = 2,
  REFLECTING = 3
};
  enum SourceFlags : int
  {
    NO_FLAGS_SET               = 0,
    APPLY_MATERIAL_SOURCE      = (1 << 0),
    APPLY_WGS_SCATTER_SOURCE   = (1 << 1),
    APPLY_AGS_SCATTER_SOURCE   = (1 << 2),
    APPLY_WGS_FISSION_SOURCE   = (1 << 3),
    APPLY_AGS_FISSION_SOURCE   = (1 << 4)
  };

  inline SourceFlags operator|(const SourceFlags f1,
                               const SourceFlags f2)
  {
    return static_cast<SourceFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
  }





//################################################################### Class def
/**A neutral particle transport solver.*/
class SteadySolver : public chi_physics::Solver
{
  typedef chi_mesh::sweep_management::CellFaceNodalMapping CellFaceNodalMapping;
protected:
  size_t source_event_tag=0;

public:
  double last_restart_write=0.0;
  lbs::Options options;

  size_t num_moments;
  size_t num_groups;
  size_t num_precursors;
  size_t max_precursors_per_material;

  std::vector<LBSGroup> groups;
  std::vector<LBSGroupset> groupsets;

  std::vector<PointSource> point_sources;

  std::map<int,std::shared_ptr<chi_physics::TransportCrossSections>> matid_to_xs_map;
  std::map<int,std::shared_ptr<chi_physics::IsotropicMultiGrpSource>> matid_to_src_map;

  std::shared_ptr<chi_math::SpatialDiscretization> discretization = nullptr;
  chi_mesh::MeshContinuumPtr grid;
  std::vector<CellFaceNodalMapping> grid_nodal_mappings;
  std::vector<UnitCellMatrices> unit_cell_matrices;
  std::vector<lbs::CellLBSView> cell_transport_views;

  //Boundaries are manipulated in chi_sweepbuffer.cc:InitializeLocalAndDownstreamBuffers
  //A default 0.0 incident boundary is loaded at the back of
  //the stack to use as default. This is loaded during initparrays
  std::vector<std::pair<BoundaryType, int>>         boundary_types;
  std::vector<std::vector<double>>                  incident_P0_mg_boundaries;
  std::vector<double>                               zero_boundary;
  std::vector<std::shared_ptr<SweepBndry>>          sweep_boundaries;

  chi_math::UnknownManager flux_moments_uk_man;

  size_t max_cell_dof_count = 0;
  uint64_t local_node_count = 0;
  uint64_t glob_node_count = 0;

  Vec phi_new = nullptr, phi_old = nullptr, q_fixed = nullptr;
  std::vector<double> q_moments_local, ext_src_moments_local;
  std::vector<double> phi_new_local, phi_old_local;
  std::vector<double> delta_phi_local;
  std::vector<std::vector<double>> psi_new_local;
  std::vector<double> precursor_new_local;

 public:
  SteadySolver (const SteadySolver&) = delete;
  SteadySolver& operator= (const SteadySolver&) = delete;

  //00
  explicit SteadySolver(const std::string& in_text_name);
  ~SteadySolver() override =default;

  //01
  void Initialize() override;
  //01a
  virtual void PerformInputChecks();
  //01b
  void PrintSimHeader();
  //01c
  void InitMaterials(std::set<int> &material_ids);
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
  //01i
  static
  bool CheckPointInsideCell(const chi_mesh::Cell& cell,
                            const chi_mesh::MeshContinuum& grid_ref,
                            const chi_mesh::Vector3& point);
  void InitializePointSources();






  //02
  void Execute() override;
  void SolveGroupset(LBSGroupset& groupset);

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
  void ExecuteWGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleWGDSADeltaPhiVector(LBSGroupset& groupset,
                                   const double *ref_phi_old,
                                   const double *ref_phi_new);
  void DisAssembleWGDSADeltaPhiVector(LBSGroupset& groupset,
                                      double *ref_phi_new);
  static void CleanUpWGDSA(LBSGroupset& groupset);

  //03e
  void InitTGDSA(LBSGroupset& groupset);
  void ExecuteTGDSA(LBSGroupset& groupset,
                    const std::vector<double>& ref_phi_old,
                    std::vector<double>& ref_phi_new);
  void AssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                   const double *ref_phi_old,
                                   const double *ref_phi_new);
  void DisAssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                      double *ref_phi_new);
  static void CleanUpTGDSA(LBSGroupset& groupset);
  //03f
  void ResetSweepOrderings(LBSGroupset& groupset);

  //04 File IO
  //04a
  void WriteRestartData(std::string folder_name, std::string file_base);
  void ReadRestartData(std::string folder_name, std::string file_base);

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

  //IterativeMethods
  virtual void SetSource(LBSGroupset& groupset,
                         std::vector<double>&  destination_q,
                         SourceFlags source_flags);
  double ComputePiecewiseChange(LBSGroupset& groupset);
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);
  bool ClassicRichardson(LBSGroupset& groupset,
                         MainSweepScheduler& sweep_scheduler,
                         SourceFlags source_flags,
                         bool log_info = true);
  bool GMRES(LBSGroupset& groupset,
             MainSweepScheduler& sweep_scheduler,
             SourceFlags lhs_src_scope,
             SourceFlags rhs_src_scope,
             bool log_info = true);

  //Vector assembly
  void SetPETScVecFromSTLvector(LBSGroupset& groupset, Vec x,
                                const std::vector<double>& y,
                                bool with_delayed_psi=false);
  void SetSTLvectorFromPETScVec(LBSGroupset& groupset, Vec x_src,
                                std::vector<double>& y,
                                bool with_delayed_psi=false);
  void ScopedCopySTLvectors(LBSGroupset& groupset,
                            const std::vector<double>& x_src,
                            std::vector<double>& y,
                            bool with_delayed_psi=false);

  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();

  //precursors
  void ComputePrecursors();
};

}

#endif
