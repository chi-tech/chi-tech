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
#include <ChiMath/chi_math.h>

#define MC_NUM_PARTICLES     1
#define MC_TFC_UPDATE_INTVL  2
#define MC_MONOENERGETIC     3
#define MC_SCATTERING_ORDER  4
#define MC_FORCE_ISOTROPIC   5
#define MC_GROUP_BOUNDS      6
#define MC_TALLY_MERGE_INTVL 7
#define MC_TALLY_MULTIPLICATION_FACTOR 8

//######################################################### Class def
/**Monte Carlo neutron particle solver.*/
class chi_montecarlon::Solver : public chi_physics::Solver
{
private:
  chi_mesh::MeshContinuum*              grid;
  SpatialDiscretization_FV*                fv_discretization;

  std::vector<int>                      matid_xs_map;

  std::vector<unsigned long long>       batch_sizes;
  std::vector<unsigned long long>       batch_sizes_per_loc;

  int                                   num_grps;

  std::vector<double>                   phi_tally_contrib;
  std::vector<double>                   phi_tally;
  std::vector<double>                   phi_tally_sqr;

  std::vector<double>                   phi_global;
  std::vector<double>                   phi_global_tally_sqr;

  std::vector<double>                   phi_local_relsigma;

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
  void ComputeTallySqr();
  void ComputeRelativeStdDev();
  void RendesvouzTallies();
};

#endif