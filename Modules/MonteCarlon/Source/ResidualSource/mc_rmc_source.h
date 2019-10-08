#ifndef _mc_rmc_source_h
#define _mc_rmc_source_h

#include "../mc_base_source.h"

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

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
  double total_Cresidual;
  double total_Lresidual;
  double total_Rresidual;

  const bool sample_uniformly;

  chi_math::CDFSampler* residual_sampler;
public:
  ResidualSource(
    chi_physics::FieldFunction* in_resid_ff,
    bool use_uniform_sampling=false);
  void Initialize(chi_mesh::MeshContinuum* ref_grid,
                  SpatialDiscretization_FV*   ref_fv_sdm);
  chi_montecarlon::Particle
  CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);

  chi_montecarlon::Particle
  UniformSampling(chi_montecarlon::RandomNumberGenerator* rng);

  chi_montecarlon::Particle
  DirectSampling(chi_montecarlon::RandomNumberGenerator* rng);
};



#endif