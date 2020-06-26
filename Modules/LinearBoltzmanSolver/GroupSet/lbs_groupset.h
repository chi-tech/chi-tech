#ifndef _lbs_groupset_h
#define _lbs_groupset_h

#include "lbs_group.h"
#include "../IterativeMethods/lbs_iterativemethods.h"

#include <ChiMath/Quadratures/LegendrePoly/legendrepoly.h>
#include <ChiMath/Quadratures/angular_quadrature_base.h>
#include <ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h>

#include <ChiPhysics/chi_physics_namespace.h>

namespace LinearBoltzman
{
  enum class AngleAggregationType
  {
    SINGLE = 1,
    POLAR = 2
  };
}

typedef chi_mesh::sweep_management::AngleAggregation AngleAgg;
typedef std::pair<int,int> GsSubSet;
typedef std::pair<int,int> AngSubSet;

#include <vector>

//################################################################### Class def
/**Group set functioning as a collection of groups*/
class LBSGroupset
{
public:
  std::vector<LBSGroup*>                       groups;
  std::shared_ptr<chi_math::AngularQuadrature> quadrature;
  std::vector<std::vector<double>>             d2m_op;
  std::vector<std::vector<double>>             m2d_op;
  chi_mesh::sweep_management::AngleAggregation* angle_agg;
  int                                          master_num_grp_subsets;
  int                                          master_num_ang_subsets;
  std::vector<GsSubSet>                        grp_subsets;
  std::vector<int>                             grp_subset_sizes;
  std::vector<AngSubSet>                       ang_subsets_top;
  std::vector<int>                             ang_subset_sizes_top;
  std::vector<AngSubSet>                       ang_subsets_bot;
  std::vector<int>                             ang_subset_sizes_bot;

  int                                          iterative_method;
  LinearBoltzman::AngleAggregationType         angleagg_method;
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
  std::vector<int>                             wgdsa_cell_dof_array_address;

  bool                                         log_sweep_events;

  double                                       latest_convergence_metric;

  //npt_groupset.cc
       LBSGroupset();
  void BuildDiscMomOperator(int scatt_order);
  void BuildMomDiscOperator(int scatt_order);
  void BuildSubsets();
public:
  void PrintSweepInfoFile(size_t ev_tag,std::string file_name);
};

#endif