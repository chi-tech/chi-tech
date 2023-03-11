#include "ss_src_context.h"

#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

namespace lbs
{

double SteadyStateSourceContext::AddSourceMoments() const
{
  return fixed_src_moments_[g_];
}

double SteadyStateSourceContext::
  AddScattering(const chi_math::SparseMatrix& S,
                const double* phi) const
{
  double value = 0.0;
  //==================== Add Across GroupSet Scattering (AGS)
  if (apply_ags_scatter_src_)
    for (const auto& [_, gp, sigma_sm] : S.Row(g_))
      if (gp < gs_i_ or gp > gs_f_)
        value += sigma_sm * phi[gp];

  //==================== Add Within GroupSet Scattering (WGS)
  if (apply_wgs_scatter_src_)
    for (const auto& [_, gp, sigma_sm] : S.Row(g_))
      if (gp >= gs_i_ and gp <= gs_f_)
      {
        if (suppress_wg_scatter_src_ and g_ == gp) continue;
          value += sigma_sm * phi[gp];
      }

  return value;
}

double SteadyStateSourceContext::
  AddPromptFission(const std::vector<double>& F_g,
                   const double* phi) const
{
  double value = 0.0;
  //==================== Add Across GroupSet Fission (AGS)
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        value += F_g[gp] * phi[gp];

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      value += F_g[gp] * phi[gp];

  return value;
}

double SteadyStateSourceContext::
  AddDelayedFission(const PrecursorList &precursors,
                    const std::vector<double> &nu_delayed_sigma_f,
                    const double *phi) const
{
  double value = 0.0;
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        for (const auto& precursor : precursors)
          value += precursor.emission_spectrum[g_] *
                   precursor.fractional_yield *
                   nu_delayed_sigma_f[gp] * phi[gp];

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
        value += precursor.emission_spectrum[g_] *
                 precursor.fractional_yield *
                 nu_delayed_sigma_f[gp] * phi[gp];

  return value;
}


}//namespace lbs