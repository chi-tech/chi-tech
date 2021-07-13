#ifndef LBS_GROUPSET_H
#define LBS_GROUPSET_H

#include "lbs_group.h"
#include "../IterativeMethods/lbs_iterativemethods.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "../lbs_structs.h"

#include "ChiPhysics/chi_physics_namespace.h"

namespace LinearBoltzmann
{
  enum class AngleAggregationType
  {
    UNDEFINED = 0,
    SINGLE = 1,
    POLAR = 2,
    AZIMUTHAL = 3,
  };
}

typedef std::pair<int,int> GsSubSet;
typedef std::pair<int,int> AngSubSet;

#include <vector>

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
  int                                          master_num_grp_subsets;
  int                                          master_num_ang_subsets;
  std::vector<GsSubSet>                        grp_subsets;
  std::vector<int>                             grp_subset_sizes;
  std::vector<AngSubSet>                       ang_subsets_top;
  std::vector<int>                             ang_subset_sizes_top;
  std::vector<AngSubSet>                       ang_subsets_bot;
  std::vector<int>                             ang_subset_sizes_bot;

  LinearBoltzmann::IterativeMethod             iterative_method;
  LinearBoltzmann::AngleAggregationType        angleagg_method;
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
  size_t                                       num_psi_unknowns_local=0;

  //lbs_groupset.cc
  explicit LBSGroupset(int in_id);
  LBSGroupset() : LBSGroupset(-1) {};
  void BuildDiscMomOperator(unsigned int scattering_order,
                            LinearBoltzmann::GeometryType geometry_type);
  void BuildMomDiscOperator(unsigned int scattering_order,
                            LinearBoltzmann::GeometryType geometry_type);
  void BuildSubsets();
public:
  void PrintSweepInfoFile(size_t ev_tag,const std::string& file_name);
};

#endif
