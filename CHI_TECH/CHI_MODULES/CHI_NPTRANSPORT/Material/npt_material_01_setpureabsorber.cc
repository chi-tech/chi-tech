#include "npt_material.h"
#include <vector>

//###################################################################
/**Sets the material as a pure absorber.*/
void NPT_MATERIAL::
SetAsPureAbsorber(int number_of_groups, int L, double sigma_t_in)
{
  //============================================= Group by group xs's
  for (int g=0;g<number_of_groups; g++)
  {
    sigma_t.push_back(sigma_t_in);
    sigma_s.push_back(0.0);
    sigma_a.push_back(sigma_t_in);
  }

  //============================================= Moment by moment transfer
  for (int ell=0; ell<=L; ell++)
  {
    NPT_MATERIAL_MOMENT_TRANFER* moment_transfer_mats =
      new NPT_MATERIAL_MOMENT_TRANFER;

    for (int g=0; g<number_of_groups; g++)
    {
      double* gprime_to_g = new double[number_of_groups];

      for (int gprime=0;gprime<number_of_groups;gprime++)
      {
        gprime_to_g[gprime] = 0.0;
      }
      moment_transfer_mats->sigma_sm_gp_to_g.push_back(gprime_to_g);
    }

    sigma_sm_matrices.push_back(moment_transfer_mats);
  }
}