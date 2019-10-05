#ifndef _lbs_linearboltzmansolver_h
#define _lbs_linearboltzmansolver_h

#include "ChiPhysics/SolverBase/chi_solver.h"

#include "GroupSet/lbs_groupset.h"
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>
#include"ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "lbs_structs.h"
#include "ChiMesh/SweepUtilities/chi_sweep.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#include <petscksp.h>

typedef chi_mesh::SweepManagement::SweepChunk SweepChunk;
typedef chi_mesh::SweepManagement::SweepScheduler MainSweepScheduler;

namespace LinearBoltzman
{
struct BoundaryTypes
{
  static const int VACUUM = 1;
  static const int INCIDENT_ISOTROPIC = 2;
};
struct SourceFlags
{
  static const bool USE_MATERIAL_SOURCE = true;
  static const bool USE_DLINV_SOURCE = false;
  static const bool SUPPRESS_PHI_OLD = true;
};




//################################################################### Class def
/**A neutral particle transport solver.*/
class Solver : public chi_physics::Solver {
 public:
  LinearBoltzman::Options options;    //In chi_npt_structs.h

  int num_moments;

  std::vector<LBSGroup*> groups;
  std::vector<LBSGroupset *> group_sets;
  std::vector<chi_physics::TransportCrossSections *> material_xs;
  std::vector<chi_physics::IsotropicMultiGrpSource *> material_srcs;
  std::vector<int> matid_to_xs_map;
  std::vector<int> matid_to_src_map;

  SpatialDiscretization *discretization;
  chi_mesh::MeshContinuum *grid;
  std::vector<LinearBoltzman::CellViewBase *> cell_transport_views;

  std::vector<int> local_cell_indices;

  //Boundaries are manipulated in chi_sweepbuffer.cc:InitializeBuffers
  //A default 0.0 incident boundary is loaded at the back of
  //the stack to use as default. This is loaded during initparrays
  std::vector<std::pair<int, int>> boundary_types;
  std::vector<std::vector<double>> incident_P0_mg_boundaries;
  std::vector<chi_mesh::SweepManagement::SPDS *> sweep_orderings;

  ChiMPICommunicatorSet comm_set;

  int max_cell_dof_count;
  unsigned long long local_dof_count;
  unsigned long long glob_dof_count;

  Vec phi_new, phi_old, q_fixed;
  std::vector<double> q_moments_local;
  std::vector<double> phi_new_local, phi_old_local, phi_oldcycle_local;
  std::vector<double> delta_phi_local;

  std::vector<int> local_cell_phi_dof_array_address;
  std::vector<int> local_cell_dof_array_address;

 public:
  //00
  Solver();
  //01
  void Initialize();
  //01a
  void ComputeNumberOfMoments();
  //01b
  void InitMaterials(std::set<int> &material_ids);
  //01c
  int InitializeParrays();
  //01d
  void InitializeCommunicators();
  //02
  void Execute();
  void SolveGroupset(int group_set_num);

  //03a
  void ComputeSweepOrderings(LBSGroupset *groupset);
  //03b
  void InitFluxDataStructures(LBSGroupset *groupset);
  //03c
  void InitAngleAggPolar(LBSGroupset *groupset);
  void InitAngleAggSingle(LBSGroupset *groupset);
  //03d
  void InitWGDSA(LBSGroupset *groupset);
  void AssembleWGDSADeltaPhiVector(LBSGroupset *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleWGDSADeltaPhiVector(LBSGroupset *groupset,
                                      double *ref_phi_new);
  //04d
  void InitTGDSA(LBSGroupset *groupset);
  void AssembleTGDSADeltaPhiVector(LBSGroupset *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleTGDSADeltaPhiVector(LBSGroupset *groupset,
                                      double *ref_phi_new);

  //04c
  void ResetSweepOrderings(LBSGroupset *groupset);

  //IterativeMethods
  void SetSource(int group_set_num,
                 bool apply_mat_src = false,
                 bool suppress_phi_old = false);
  double ComputePiecewiseChange(LBSGroupset *groupset);
  SweepChunk *SetSweepChunk(int group_set_num);
  void ClassicRichardson(int group_set_num);
  void GMRES(int group_set_num);
  void AssembleVectors(LBSGroupset *groupset);
  void AssembleVector(LBSGroupset *groupset, Vec x, double *y);
  void DisAssembleVector(LBSGroupset *groupset, Vec x_src, double *y);
  void DisAssembleVectorLocalToLocal(LBSGroupset *groupset, double *x_src, double *y);
  void ConvergeCycles(MainSweepScheduler& sweepScheduler,
                      SweepChunk* sweep_chunk,
                      LBSGroupset *groupset);
};

}

#endif