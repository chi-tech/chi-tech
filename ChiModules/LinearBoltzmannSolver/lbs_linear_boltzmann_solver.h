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

#include <petscksp.h>

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

namespace LinearBoltzmann
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
class Solver : public chi_physics::Solver
{
  typedef chi_mesh::sweep_management::CellFaceNodalMapping CellFaceNodalMapping;
protected:
  size_t source_event_tag=0;

public:
  double last_restart_write=0.0;
  LinearBoltzmann::Options options;    //In chi_npt_structs.h

  size_t num_moments;
  size_t num_groups;
  size_t num_precursors;
  size_t max_precursors_per_material;

  std::vector<LBSGroup> groups;
  std::vector<LBSGroupset> groupsets;
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> material_xs;
  std::vector<std::shared_ptr<chi_physics::IsotropicMultiGrpSource>> material_srcs;
  std::vector<int> matid_to_xs_map;
  std::vector<int> matid_to_src_map;

  std::shared_ptr<SpatialDiscretization> discretization = nullptr;
  chi_mesh::MeshContinuumPtr grid;
  std::vector<CellFaceNodalMapping> grid_nodal_mappings;
  std::vector<LinearBoltzmann::CellLBSView> cell_transport_views;

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
  //00
  explicit Solver(const std::string& in_text_name);
  ~Solver() override =default;

  //01
  void Initialize() override;
  //01a
  virtual void PerformInputChecks();
  void ComputeNumberOfMoments();
  void PrintSimHeader();
  //01b
  void InitMaterials(std::set<int> &material_ids);
  //01c
  virtual void InitializeSpatialDiscretization();
  //01c
  void InitializeBoundaries();
  //01d
  virtual void InitializeParrays();
  //01e
  void InitializeGroupsets();

  //02
  void Execute() override;
  void SolveGroupset(LBSGroupset& groupset);

  //03a
  void ComputeSweepOrderings(LBSGroupset& groupset) const;
  void ComputeSweepOrderingsAngleAggSingle(LBSGroupset& groupset) const;
  void ComputeSweepOrderingsAngleAggPolar(LBSGroupset& groupset) const;
  void ComputeSweepOrderingsAngleAggAzimuthal(LBSGroupset& groupset) const;
  //03b
  void InitFluxDataStructures(LBSGroupset& groupset);
  //03c
  void InitAngleAggPolar(LBSGroupset& groupset);
  void InitAngleAggSingle(LBSGroupset& groupset);
  void InitAngleAggAzimuthal(LBSGroupset& groupset);
  //03d
  void InitWGDSA(LBSGroupset& groupset);
  void AssembleWGDSADeltaPhiVector(LBSGroupset& groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleWGDSADeltaPhiVector(LBSGroupset& groupset,
                                      double *ref_phi_new);
  void CleanUpWGDSA(LBSGroupset& groupset);
  //03e
  void InitTGDSA(LBSGroupset& groupset);
  void AssembleTGDSADeltaPhiVector(LBSGroupset& groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                      double *ref_phi_new);
  void CleanUpTGDSA(LBSGroupset& groupset);
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
                                bool with_delayed_psi= false);
  void SetSTLvectorFromPETScVec(LBSGroupset& groupset, Vec x_src,
                                std::vector<double>& y,
                                bool with_delayed_psi= false);
  void ScopedCopySTLvectors(LBSGroupset& groupset,
                            const std::vector<double>& x_src,
                            std::vector<double>& y);

  //compute_balance
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);
  void ComputeBalance();

  //precursors
  void ComputePrecursors();
};

}

#endif
