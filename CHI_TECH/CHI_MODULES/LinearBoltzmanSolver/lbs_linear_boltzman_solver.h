#ifndef _lbs_linearboltzmansolver_h
#define _lbs_linearboltzmansolver_h

#include "../../ChiPhysics/SolverBase/chi_solver.h"

#include "CHI_MODULES/LinearBoltzmanSolver/GroupSet/lbs_groupset.h"
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>
#include"ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "lbs_structs.h"
#include "../../CHI_MESH/CHI_SWEEP/chi_sweep.h"
#include "../../ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#include <petscksp.h>

typedef chi_mesh::SweepManagement::SweepChunk SweepChunk;

#define VACUUM             301
#define INCIDENT_ISOTROPIC 302

#define USE_MATERIAL_SOURCE true
#define USE_DLINV_SOURCE false
#define SUPPRESS_PHI_OLD true

//################################################################### Class def
/**A neutral particle transport solver.*/
class LinearBoltzmanSolver : public chi_physics::Solver
{
public:
  LBS_OPTIONS options;    //In chi_npt_structs.h

  int num_moments;

  std::vector<LBS_GROUP*>                            groups;
  std::vector<LBS_GROUPSET*>                         group_sets;
  std::vector<chi_physics::TransportCrossSections*>  material_xs;
  std::vector<chi_physics::IsotropicMultiGrpSource*> material_srcs;
  std::vector<int>                                   matid_to_xs_map;
  std::vector<int>                                   matid_to_src_map;


  SpatialDiscretization*                                discretization;
  chi_mesh::MeshContinuum*                           grid;
  std::vector<LBS_CELLVIEW*>                         cell_transport_views;

  std::vector<int>                                   local_cell_indices;

  //Boundaries are manipulated in chi_sweepbuffer.cc:InitializeBuffers
  //A default 0.0 incident boundary is loaded at the back of
  //the stack to use as default. This is loaded during initparrays
  std::vector<std::pair<int,int>>                    boundary_types;
  std::vector<std::vector<double>>                   incident_P0_mg_boundaries;
  std::vector<chi_mesh::SweepManagement::SPDS*>      sweep_orderings;

  CHI_MPI_COMMUNICATOR_SET                           comm_set;

  int max_cell_dof_count;
  unsigned long long local_dof_count;
  unsigned long long glob_dof_count;

  Vec phi_new, phi_old, q_fixed;
  std::vector<double> q_fixed_local;
  std::vector<double> q_moments_local;
  std::vector<double> phi_new_local, phi_old_local;
  std::vector<double> delta_phi_local;

  std::vector<int> local_indices;
  std::vector<int> global_indices;

  std::vector<int> local_cell_phi_dof_array_address;
  std::vector<int> local_cell_dof_array_address;

  ISLocalToGlobalMapping ltog;

public:
  //00
       LinearBoltzmanSolver();
  //01
  void Initialize();
  //01a
  void SetPartitioning();
  void ComputeNumberOfMoments();
  //01b
  void InitMaterials(std::set<int>& material_ids);
  //01c
  int InitializeParrays();
  //01d
  void InitializeCommunicators();
  //02
  void Execute();
  void SolveGroupset(int group_set_num);

  //03a
  void ComputeSweepOrderings(LBS_GROUPSET *groupset);
  //03b
  void InitFluxDataStructures(LBS_GROUPSET *groupset);
  //03c
  void InitAngleAggPolar(LBS_GROUPSET *groupset);
  void InitAngleAggSingle(LBS_GROUPSET *groupset);
  //03d
  void InitWGDSA(LBS_GROUPSET *groupset);
  void AssembleWGDSADeltaPhiVector(LBS_GROUPSET *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleWGDSADeltaPhiVector(LBS_GROUPSET *groupset,
                                      double *ref_phi_new);
  //04d
  void InitTGDSA(LBS_GROUPSET* groupset);
  void AssembleTGDSADeltaPhiVector(LBS_GROUPSET *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleTGDSADeltaPhiVector(LBS_GROUPSET *groupset,
                                      double *ref_phi_new);


  //04c
  void ResetSweepOrderings(LBS_GROUPSET *groupset);



  //IterativeMethods
  void        SetSource(int group_set_num,
                        bool apply_mat_src=false,
                        bool suppress_phi_old=false);
  double      ComputePiecewiseChange(LBS_GROUPSET* groupset);
  SweepChunk* SetSweepChunk(int group_set_num);
  void        ClassicRichardson(int group_set_num);
  void        GMRES(int group_set_num);
  void        AssembleVectors(LBS_GROUPSET *groupset);
  void        AssembleVector(LBS_GROUPSET *groupset, Vec x, double *y);
  void        DisAssembleVector(LBS_GROUPSET *groupset, Vec x_src, double *y);
  void        DisAssembleVectorLocalToLocal(LBS_GROUPSET *groupset, double* x_src, double *y);
};

#endif