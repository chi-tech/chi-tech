#ifndef LBS_GROUPSET_H
#define LBS_GROUPSET_H

#include "lbs_group.h"
#include "LBSSteadyState/IterativeMethods/lbs_iterativemethods.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "LBSSteadyState/lbs_structs.h"
#include "LBSSteadyState/lbs_make_subset.h"

#include "ChiPhysics/chi_physics_namespace.h"

#include "LBSSteadyState/Acceleration/acceleration.h"

namespace lbs::acceleration
{
  class DiffusionMIPSolver;
}

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
  std::shared_ptr<chi_math::AngularQuadrature> quadrature = nullptr;
  chi_mesh::sweep_management::AngleAggregation angle_agg;
  std::vector<SPDS_ptr>                        sweep_orderings;
  UniqueSOGroupings                            unique_so_groupings;
  DirIDToSOMap                                 dir_id_to_so_map;

  int                                          master_num_grp_subsets = 1;
  int                                          master_num_ang_subsets = 1;

  std::vector<SubSetInfo>                      grp_subset_infos;

  IterativeMethod      iterative_method = IterativeMethod::CLASSICRICHARDSON;
  AngleAggregationType angleagg_method = AngleAggregationType::POLAR;
  double               residual_tolerance = 1.0e-6;
  int                  max_iterations = 200;
  int                  gmres_restart_intvl = 30;

  bool                 allow_cycles = false;
  bool                 log_sweep_events = false;

  bool                 apply_wgdsa = false;
  bool                 apply_tgdsa = false;
  int                  wgdsa_max_iters = 30;
  int                  tgdsa_max_iters = 30;
  double               wgdsa_tol = 1.0e-4;
  double               tgdsa_tol = 1.0e-4;
  bool                 wgdsa_verbose = false;
  bool                 tgdsa_verbose = false;
  std::string          wgdsa_string;
  std::string          tgdsa_string;

  std::shared_ptr<lbs::acceleration::DiffusionMIPSolver> wgdsa_solver;
  std::shared_ptr<lbs::acceleration::DiffusionMIPSolver> tgdsa_solver;

  struct TwoGridAccelerationInfo
  {
    std::map<int, acceleration::TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    acceleration::EnergyCollapseScheme scheme =
      acceleration::EnergyCollapseScheme::JFULL;
  }tg_acceleration_info;


  chi_math::UnknownManager psi_uk_man;

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
