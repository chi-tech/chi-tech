#include "mc_base_source.h"
#include<math.h>

#include "../RandomNumberGenerator/montecarlon_rng.h"

//#########################################################
/**Default constructor*/
chi_montecarlon::Source::Source()
{
  type_index = MC_BASE_SRC;
}

void chi_montecarlon::Source::Initialize(chi_mesh::MeshContinuum* ref_grid,
                                         SpatialDiscretization_FV*   ref_fv_sdm)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
}

//#########################################################
/**Create a default particle from a point.*/
chi_montecarlon::Particle chi_montecarlon::Source::CreateParticle(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  new_particle.pos.x = 0.0;
  new_particle.pos.y = 0.0;
  new_particle.pos.z = 0.0;

  double costheta = rng->Rand()*2.0 - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  new_particle.dir.x = sin(theta)*cos(varphi);
  new_particle.dir.y = sin(theta)*sin(varphi);
  new_particle.dir.z = cos(theta);

  new_particle.egrp = 0;
  new_particle.w = 1.0;


  return new_particle;
}