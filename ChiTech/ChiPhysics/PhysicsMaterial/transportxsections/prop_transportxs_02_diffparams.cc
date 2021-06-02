#include "material_property_transportxsections.h"

#include <ChiLog/chi_log.h>

extern ChiLog& chi_log;

void chi_physics::TransportCrossSections::ComputeDiffusionParameters()
{
  if (diffusion_initialized)
    return;

  diffg.resize(num_groups, 1.0);
  sigma_s_gtog.resize(num_groups, 0.0);
  sigma_rg.resize(num_groups, 0.1);
  for (int g=0; g < num_groups; g++)
  {
    //====================================== Determine transport correction
    double sigs_g_1 = 0.0;
    if (transfer_matrix.size()>1)
    {
      for (int gp=0; gp < num_groups; gp++)
      {
        size_t num_cols = transfer_matrix[1].rowI_indices[gp].size();
        for (int j=0; j<num_cols; j++)
        {
          if (transfer_matrix[1].rowI_indices[gp][j] == g)
          {
            sigs_g_1 += transfer_matrix[1].rowI_values[gp][j];
            break;
          }
        }//for j
      }//for gp
    }//if moment 1 available

    //====================================== Determine diffcoeff
    if (sigs_g_1 >= sigma_tg[g])
    {
      sigs_g_1 = 0.0;
      chi_log.Log(LOG_0WARNING)
        << "Transport corrected diffusion coefficient failed for group "
        << g << " in call to "
        << "chi_physics::TransportCrossSections::ComputeDiffusionParameters."
        << " sigma_t=" << sigma_tg[g] << " sigs_g_(m=1)=" << sigs_g_1;
    }
    diffg[g] = (fmin(1.0e12,1.0/3.0/(sigma_tg[g]-sigs_g_1)));
    //diffg[g] = 1.0/3.0/(sigma_tg[g]-sigs_g_1);

    //====================================== Determine in group scattering
    size_t num_cols = transfer_matrix[0].rowI_indices[g].size();
    for (int j=0; j<num_cols; j++)
    {
      if (transfer_matrix[0].rowI_indices[g][j] == g)
      {
        sigma_s_gtog[g] = transfer_matrix[0].rowI_values[g][j];
        break;
      }
    }

    //====================================== Determine removal cross-section
    sigma_rg[g] = std::max(0.0,sigma_tg[g] - sigma_s_gtog[g]);
  }//for g

  //======================================== Compute two grid energy collapse
  chi_log.Log(LOG_0) << "Performing Energy collapse.";
  EnergyCollapse(xi_Jfull_g,D_jfull,sigma_a_jfull,E_COLLAPSE_JACOBI);
  EnergyCollapse(xi_Jpart_g,D_jpart,sigma_a_jpart,E_COLLAPSE_PARTIAL_JACOBI);

  diffusion_initialized = true;
}