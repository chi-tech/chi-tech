#ifndef _solver_montecarlon
#define _solver_montecarlon
#include<iostream>

#include"../chi_montecarlon.h"
#include"ChiPhysics/SolverBase/chi_solver.h"
#include "../RandomNumberGenerator/montecarlon_rng.h"
#include "../Source/mc_base_source.h"
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <FiniteVolume/fv.h>
#include <PiecewiseLinear/pwl.h>
#include <ChiMath/chi_math.h>

namespace chi_montecarlon
{
  enum Property{
    NUM_PARTICLES     = 1,
    TFC_UPDATE_INTVL  = 2,
    MONOENERGETIC     = 3,
    SCATTERING_ORDER  = 4,
    FORCE_ISOTROPIC   = 5,
    GROUP_BOUNDS      = 6,
    TALLY_MERGE_INTVL = 7,
    TALLY_MULTIPLICATION_FACTOR = 8,
    MAKE_PWLD_SOLUTION = 9
  };
}

//######################################################### Class def
/**Monte Carlo neutron particle solver.*/
class chi_montecarlon::Solver : public chi_physics::Solver
{
private:
  chi_mesh::MeshContinuum*              grid;
  SpatialDiscretization_FV*             fv_discretization;
  SpatialDiscretization_PWL*            pwl_discretization;

  std::vector<int>                      matid_xs_map;

  std::vector<unsigned long long>       batch_sizes;
  std::vector<unsigned long long>       batch_sizes_per_loc;

  int                                   num_grps;

  //FV tallies
  std::vector<double>                   phi_tally_contrib;
  std::vector<double>                   phi_tally;
  std::vector<double>                   phi_tally_sqr;

  std::vector<double>                   phi_global;
  std::vector<double>                   phi_global_tally_sqr;

  std::vector<double>                   phi_local_relsigma;

  //PWL tallies
  int                                   num_moms;
  std::vector<double>                   phi_pwl_tally_contrib;
  std::vector<double>                   phi_pwl_tally;
  std::vector<double>                   phi_pwl_tally_sqr;

  std::vector<double>                   phi_pwl_global;
  std::vector<double>                   phi_pwl_global_tally_sqr;

  std::vector<double>                   phi_pwl_local_relsigma;

  std::vector<int>                      local_cell_pwl_dof_array_address;


  //runtime quantities
  int                                   current_batch;
  unsigned long long                    nps;
  unsigned long long                    nps_global;
  double                                max_relative_error;
  double                                max_relative_error2;
  double                                max_relative_error3;
  chi_math::CDFSampler*                 surface_source_sampler;


public:
  RandomNumberGenerator                 rng0;
  std::vector<chi_montecarlon::Source*> sources;

  //options
  double                                tolerance;
  unsigned long long                    num_particles;
  int                                   tfc_update_interval;
  bool                                  mono_energy;
  int                                   scattering_order;
  bool                                  force_isotropic;
  int                                   group_hi_bound;
  int                                   group_lo_bound;
  unsigned long long                    tally_rendezvous_intvl;
  double                                tally_multipl_factor;
  bool                                  make_pwld;


  //derived from options-set during init


private:
  chi_mesh::EmptyRegion empty_region;

public:
  //00
       Solver();
  //01
  bool Initialize();

  //02
  void Execute();

private:
  //03
  void Raytrace(chi_montecarlon::Particle* prtcl);

  //04
  std::pair<int,chi_mesh::Vector>
  ProcessScattering(chi_montecarlon::Particle* prtcl,
                    chi_physics::TransportCrossSections* xs);
  //05
  void ContributeTally(chi_montecarlon::Particle* prtcl,
                       chi_mesh::Vector pf);
  void RendesvouzTallies();
  void RendesvouzPWLTallies();

  void ComputeRelativeStdDev();

  void NormalizeTallies();
  void NormalizePWLTallies();




};

#endif