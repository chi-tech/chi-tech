#ifndef _chi_nptransport_h
#define _chi_nptransport_h

#include "../../CHI_PHYSICS/CHI_SOLVER/chi_solver.h"

#include "GroupSet/npt_groupset.h"
#include <CHI_PHYSICS/CHI_PHYSICSMATERIAL/property10_transportxsections.h>
#include <CHI_PHYSICS/CHI_PHYSICSMATERIAL/property11_isotropic_mg_src.h>
#include"../../CHI_MATH/CHI_DISCRETIZATION/chi_discretization.h"
#include "chi_npt_structs.h"
#include "../../CHI_MESH/CHI_SWEEP/chi_sweep.h"
#include "../../CHI_MATH/SparseMatrix/chi_math_sparse_matrix.h"

#include <petscksp.h>

typedef chi_mesh::SweepManagement::SweepChunk SweepChunk;

#define VACUUM             301
#define INCIDENT_ISOTROPIC 302

#define USE_MATERIAL_SOURCE true
#define USE_DLINV_SOURCE false
#define SUPPRESS_PHI_OLD true

//################################################################### Class def
/**A neutral particle transport solver.*/
class CHI_NPTRANSPORT : public chi_physics::Solver
{
public:
  NPT_OPTIONS options;    //In chi_npt_structs.h

  int num_moments;

  std::vector<NPT_GROUP*>                            groups;
  std::vector<NPT_GROUPSET*>                         group_sets;
  std::vector<chi_physics::TransportCrossSections*>  material_xs;
  std::vector<chi_physics::IsotropicMultiGrpSource*> material_srcs;
  std::vector<int>                                   matid_to_xs_map;
  std::vector<int>                                   matid_to_src_map;


  CHI_DISCRETIZATION*                                discretization;
  chi_mesh::MeshContinuum*                           grid;
  std::vector<NPT_CELLVIEW*>                         cell_transport_views;

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
       CHI_NPTRANSPORT();
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
  void ComputeSweepOrderings(NPT_GROUPSET *groupset);
  //03b
  void InitFluxDataStructures(NPT_GROUPSET *groupset);
  //03c
  void InitAngleAggPolar(NPT_GROUPSET *groupset);
  void InitAngleAggSingle(NPT_GROUPSET *groupset);
  //03d
  void InitWGDSA(NPT_GROUPSET *groupset);
  void AssembleWGDSADeltaPhiVector(NPT_GROUPSET *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleWGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                      double *ref_phi_new);
  //04d
  void InitTGDSA(NPT_GROUPSET* groupset);
  void AssembleTGDSADeltaPhiVector(NPT_GROUPSET *groupset, double *ref_phi_old,
                                   double *ref_phi_new);
  void DisAssembleTGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                      double *ref_phi_new);


  //04c
  void ResetSweepOrderings(NPT_GROUPSET *groupset);



  //IterativeMethods
  void        SetSource(int group_set_num,
                        bool apply_mat_src=false,
                        bool suppress_phi_old=false);
  double      ComputePiecewiseChange(NPT_GROUPSET* groupset);
  SweepChunk* SetSweepChunk(int group_set_num);
  void        ClassicRichardson(int group_set_num);
  void        GMRES(int group_set_num);
  void        AssembleVectors(NPT_GROUPSET *groupset);
  void        AssembleVector(NPT_GROUPSET *groupset, Vec x, double *y);
  void        DisAssembleVector(NPT_GROUPSET *groupset, Vec x_src, double *y);
  void        DisAssembleVectorLocalToLocal(NPT_GROUPSET *groupset, double* x_src, double *y);
};

#endif