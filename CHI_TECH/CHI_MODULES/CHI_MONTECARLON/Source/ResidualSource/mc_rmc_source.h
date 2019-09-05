#ifndef _mc_rmc_source_h
#define _mc_rmc_source_h

#include "../mc_base_source.h"

#include <CHI_PHYSICS/chi_physics.h>

#include <CHI_MATH/chi_math.h>

//###################################################################
/**Residual source class.*/
class chi_montecarlon::ResidualSource : public chi_montecarlon::Source
{
private:
  chi_physics::FieldFunction* resid_ff;
  std::vector<double> cell_residuals;
  std::vector<double> cell_interior_residual;
  std::vector<double> cell_surface_residualL;
  std::vector<double> cell_surface_residualR;
  std::vector<double> cell_residual_cdf;

  double total_residual;

  chi_math::CDFSampler* residual_sampler;
public:
  ResidualSource(chi_physics::FieldFunction* in_resid_ff);
  void Initialize(chi_mesh::MeshContinuum* ref_grid,
                  CHI_DISCRETIZATION_FV*   ref_fv_sdm);
  chi_montecarlon::Particle
  CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);
};



#endif