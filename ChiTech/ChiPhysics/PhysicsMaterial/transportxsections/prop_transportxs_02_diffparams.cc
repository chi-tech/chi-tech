#include "material_property_transportxsections.h"
#include "chi_runtime.h"
#include "chi_log.h"


//######################################################################
void chi_physics::TransportCrossSections::
ComputeDiffusionParameters()
{
  if (diffusion_initialized)
    return;

  //initialize diffusion data
  diffusion_coeff.resize(num_groups, 1.0);
  sigma_s_gtog.resize(num_groups, 0.0);
  sigma_removal.resize(num_groups, 0.1);

  //perfom computations group-wise
  const auto& S = transfer_matrices;
  for (unsigned int g = 0; g < num_groups; ++g)
  {
    //============================================================
    // Determine transport correction
    //============================================================

    double sig_1 = 0.0;
    if (S.size() > 1)
    {
      for (unsigned int gp = 0; gp < num_groups; ++gp)
      {
        const auto& cols = S[1].rowI_indices[gp];
        const auto& vals = S[1].rowI_values[gp];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sig_1 += vals[t];
            break;
          }
      }//for gp
    }//if moment 1 available

    //============================================================
    // Compute diffusion coefficient
    //============================================================

    if (sig_1 >= sigma_t[g])
    {
      sig_1 = 0.0;
      chi::log.LogAllWarning()
          << "Transport corrected diffusion coefficient failed for group "
          << g << " in call to " << __FUNCTION__ << ". "
          << "sigma_t=" << sigma_t[g] << " sigs_g_(m=1)=" << sig_1
          << ". Setting sigs_g_(m=1) to zero for this group.";
    }

    //compute the diffusion coefficient
    //cap the value for when sig_t - sig_1 is near zero
    diffusion_coeff[g] = std::fmin(1.0e12,
                                   1.0 / 3.0 / (sigma_t[g] - sig_1));

    //============================================================
    // Determine within group scattering
    //============================================================

    if (!S.empty())
    {
      const auto& cols = S[0].rowI_indices[g];
      const auto& vals = S[0].rowI_values[g];
      for (size_t t = 0; t < cols.size(); ++t)
        if (cols[t] == g)
        {
          sigma_s_gtog[g] = vals[t];
          break;
        }
    }

    //============================================================
    // Compute removal cross-section
    //============================================================

    sigma_removal[g] = std::max(0.0, sigma_t[g] - sigma_s_gtog[g]);
  }//for g

  diffusion_initialized = true;
}