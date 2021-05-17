#ifndef LINEAR_BOLTZMANN_SOLVER_H
#define LINEAR_BOLTZMANN_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "GroupSet/lbs_groupset.h"
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
struct SourceFlags
{
  static const bool USE_MATERIAL_SOURCE = true;
  static const bool NO_MATERIAL_SOURCE = false;
  static const bool USE_PHI_SCATTER_SOURCE = true;
  static const bool NO_PHI_SCATTER_SOURCE = false;
};





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

  std::vector<LBSGroup> groups;
  std::vector<LBSGroupset> group_sets;
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> material_xs;
  std::vector<std::shared_ptr<chi_physics::IsotropicMultiGrpSource>> material_srcs;
  std::vector<int> matid_to_xs_map;
  std::vector<int> matid_to_src_map;

  std::shared_ptr<SpatialDiscretization> discretization;
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

  size_t max_cell_node_count;
  size_t local_node_count;
  size_t globl_node_count;

  Vec phi_new, phi_old, q_fixed;
  std::vector<double> q_moments_local;
  std::vector<double> phi_new_local, phi_old_local;
  std::vector<double> delta_phi_local;

 public:
  //00
  Solver();
  ~Solver() override =default;
  //01
  virtual void Initialize();
  //01a
  void PerformInputChecks();
  void ComputeNumberOfMoments();
  void PrintSimHeader();
  //01b
  void InitMaterials(std::set<int> &material_ids);
  //01c
  void InitializeSpatialDiscretization();
  //01c
  void InitializeBoundaries();
  //01d
  virtual void InitializeParrays();
  //01e
  void InitializeGroupsets();
  //02
  void Execute() override;
  void SolveGroupset(LBSGroupset& groupset,
                     int group_set_num);
  //02b
  void ComputeBalance();
  void SweepAllGroupSets();

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

  //04
  void WriteRestartData(std::string folder_name, std::string file_base);
  void ReadRestartData(std::string folder_name, std::string file_base);
  void WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                                  const std::string& file_base);
  void ReadGroupsetAngularFluxes(LBSGroupset& groupset,
                                 const std::string& file_base);

  //IterativeMethods
  virtual void SetSource(LBSGroupset& groupset,
                 bool apply_mat_src,
                 bool apply_phi_old_scatter_src);
  double ComputePiecewiseChange(LBSGroupset& groupset);
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);
  bool ClassicRichardson(LBSGroupset& groupset,
                         int group_set_num,
                         SweepChunk& sweep_chunk,
                         MainSweepScheduler & sweepScheduler,
                         bool log_info = true);
  bool GMRES(LBSGroupset& groupset,
             int group_set_num,
             SweepChunk& sweep_chunk,
             MainSweepScheduler & sweepScheduler,
             bool log_info = true);

  //Vector assembly
  void AssembleVector(LBSGroupset& groupset, Vec x, double *y,bool with_delayed_psi=false);
  void DisAssembleVector(LBSGroupset& groupset, Vec x_src, double *y,bool with_delayed_psi=false);
  void DisAssembleVectorLocalToLocal(LBSGroupset& groupset, double *x_src, double *y);

};

}

#endif
