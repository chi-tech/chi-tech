#ifndef LBS_GROUPSET_H
#define LBS_GROUPSET_H

#include "lbs_group.h"
#include "../IterativeMethods/lbs_iterativemethods.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "../lbs_structs.h"
#include "../lbs_make_subset.h"

#include "ChiPhysics/chi_physics_namespace.h"

namespace lbs
{

//################################################################### Class def
/**Group set functioning as a collection of groups*/
class LBSGroupset
{
protected:
  typedef std::shared_ptr<chi_mesh::sweep_management::SPDS> SPDS_ptr;
public:
  int                                          id;
  std::vector<LBSGroup>                        groups;
  std::shared_ptr<chi_math::AngularQuadrature> quadrature;
  chi_mesh::sweep_management::AngleAggregation angle_agg;
  std::vector<SPDS_ptr>                        sweep_orderings;
  UniqueSOGroupings                            unique_so_groupings;
  DirIDToSOMap                                 dir_id_to_so_map;

  int                                          master_num_grp_subsets;
  int                                          master_num_ang_subsets;

  std::vector<SubSetInfo>                      grp_subset_infos;

  IterativeMethod                              iterative_method;
  AngleAggregationType                         angleagg_method;
  double                                       residual_tolerance;
  int                                          max_iterations;
  int                                          gmres_restart_intvl;
  bool                                         apply_wgdsa;
  bool                                         apply_tgdsa;
  int                                          wgdsa_max_iters;
  int                                          tgdsa_max_iters;
  double                                       wgdsa_tol;
  double                                       tgdsa_tol;
  bool                                         wgdsa_verbose;
  bool                                         tgdsa_verbose;
  std::string                                  wgdsa_string;
  std::string                                  tgdsa_string;

  bool                                         allow_cycles;

  chi_physics::Solver*                         wgdsa_solver;
  chi_physics::Solver*                         tgdsa_solver;

  bool                                         log_sweep_events;

  chi_math::UnknownManager                     psi_uk_man;

  //lbs_groupset.cc
  explicit LBSGroupset(int in_id);
  LBSGroupset() : LBSGroupset(-1) {};
  void BuildDiscMomOperator(unsigned int scattering_order,
                            GeometryType geometry_type);
  void BuildMomDiscOperator(unsigned int scattering_order,
                            GeometryType geometry_type);
  void BuildSubsets();
public:
  void PrintSweepInfoFile(size_t ev_tag,const std::string& file_name);
};
}




#endif
