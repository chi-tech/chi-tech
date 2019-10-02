#ifndef _lbs_groupset_h
#define _lbs_groupset_h

#include "lbs_group.h"
#include "../IterativeMethods/lbs_iterativemethods.h"

#include <ChiMath/Quadratures/LegendrePoly/legendrepoly.h>
#include <ChiMath/Quadratures/product_quadrature.h>
#include <ChiMesh/SweepUtilities/chi_angleaggregation.h>

#include <ChiPhysics/chi_physics_namespace.h>

typedef chi_mesh::SweepManagement::AngleAggregation AngleAgg;

#define NPT_ANGAGG_SINGLE 1
#define NPT_ANGAGG_POLAR 2

typedef std::pair<int,int> GsSubSet;
typedef std::pair<int,int> AngSubSet;

#include <vector>

//################################################################### Class def
/**Group set functioning as a collection of groups*/
class LBSGroupset
{
public:
  std::vector<LBSGroup*>                       groups;
  chi_math::ProductQuadrature*                 quadrature;
  std::vector<std::vector<double>>             d2m_op;
  std::vector<std::vector<double>>             m2d_op;
  chi_mesh::SweepManagement::AngleAggregation* angle_agg;
  int                                          master_num_grp_subsets;
  int                                          master_num_ang_subsets;
  std::vector<GsSubSet>                        grp_subsets;
  std::vector<int>                             grp_subset_sizes;
  std::vector<AngSubSet>                       ang_subsets_top;
  std::vector<int>                             ang_subset_sizes_top;
  std::vector<AngSubSet>                       ang_subsets_bot;
  std::vector<int>                             ang_subset_sizes_bot;

  int                                          iterative_method;
  int                                          angleagg_method;
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

  //npt_groupset.cc
       LBSGroupset();
  void BuildDiscMomOperator(int scatt_order);
  void BuildMomDiscOperator(int scatt_order);
  void BuildSubsets();


};

#endif