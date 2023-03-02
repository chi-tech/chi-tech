#include "multigroup_xs.h"
#include "chi_runtime.h"
#include "chi_log.h"


//######################################################################
void chi_physics::MultiGroupXS::ComputeDiffusionParameters()
{
  if (diffusion_initialized_)
    return;

  //initialize diffusion data
  diffusion_coeff_.resize(num_groups_, 1.0);
  sigma_s_gtog_.resize(num_groups_, 0.0);
  sigma_removal_.resize(num_groups_, 0.1);

  //perfom computations group-wise
  const auto& S = transfer_matrices_;
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    //============================================================
    // Determine transport correction
    //============================================================

    double sig_1 = 0.0;
    if (S.size() > 1)
    {
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        const auto& cols = S[1].rowI_indices_[gp];
        const auto& vals = S[1].rowI_values_[gp];
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

    if (sig_1 >= sigma_t_[g])
    {
      sig_1 = 0.0;
      chi::log.Log0Warning()
          << "Transport corrected diffusion coefficient failed for group "
          << g << " in call to " << __FUNCTION__ << ". "
          << "sigma_t=" << sigma_t_[g] << " sigs_g_(m=1)=" << sig_1
          << ". Setting sigs_g_(m=1) to zero for this group.";
    }

    //compute the diffusion coefficient
    //cap the value for when sig_t - sig_1 is near zero
    diffusion_coeff_[g] = std::fmin(1.0e12,
                                   1.0 / 3.0 / (sigma_t_[g] - sig_1));

    //============================================================
    // Determine within group scattering
    //============================================================

    if (!S.empty())
    {
      const auto& cols = S[0].rowI_indices_[g];
      const auto& vals = S[0].rowI_values_[g];
      for (size_t t = 0; t < cols.size(); ++t)
        if (cols[t] == g)
        {
          sigma_s_gtog_[g] = vals[t];
          break;
        }
    }

    //============================================================
    // Compute removal cross section
    //============================================================

    sigma_removal_[g] = std::max(0.0, sigma_t_[g] - sigma_s_gtog_[g]);
  }//for g

  diffusion_initialized_ = true;
}