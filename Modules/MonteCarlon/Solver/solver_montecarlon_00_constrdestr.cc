#include "solver_montecarlon.h"

//###################################################################
/**Default constructor*/
chi_montecarlon::Solver::Solver()
{
  num_grps = 1;

  this->tolerance = 0.000001;

  //options
  num_particles = 1000;
  tfc_update_interval = 2000;
  mono_energy = false;
  scattering_order = 10;
  force_isotropic = false;
  group_hi_bound = -1;
  group_lo_bound = -1;
  tally_rendezvous_intvl = 100000;
  tally_multipl_factor = 1.0;

  max_relative_error = 0.0;

}