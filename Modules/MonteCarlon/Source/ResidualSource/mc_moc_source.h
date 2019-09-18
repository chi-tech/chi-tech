#ifndef _mc_moc_source_h
#define _mc_moc_source_h

#include "../mc_base_source.h"

#include <ChiMath/Quadratures/quadrature_gausslegendre.h>

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

//###################################################################
/**Residual source class.*/
class chi_montecarlon::ResidualMOCSource : public chi_montecarlon::Source
{
private:
  chi_physics::FieldFunction* resid_ff;

  std::vector<std::vector<double>> cell_dof_phi;
  chi_math::QuadratureGaussLegendre quadrature;

  std::vector<std::vector<double>> cell_z_i_star;
  std::vector<std::vector<double>> cell_phi_star;

  double                           total_uncollided_flux;
  std::vector<double>              cell_total_subfluxes;
  std::vector<std::vector<double>> cell_subfluxes;
  std::vector<double>              cell_cdf;
  std::vector<double>              cell_sigma_s;
  std::vector<double>              cell_sigma_t;

  const bool sample_uniformly;
  size_t num_subdivs;

  chi_math::CDFSampler* cell_sampler;
public:
  ResidualMOCSource(
    chi_physics::FieldFunction* in_resid_ff,
    bool use_uniform_sampling=false);
  void Initialize(chi_mesh::MeshContinuum* ref_grid,
                  SpatialDiscretization_FV*   ref_fv_sdm);
  chi_montecarlon::Particle
  CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);
};



#endif