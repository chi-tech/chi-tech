#include "CHI_PHYSICS/CHI_PHYSICSMATERIAL/property10_transportxsections.h"

//###################################################################
/**Default constructor.*/
chi_physics::TransportCrossSections::TransportCrossSections()
{
  type_index = TRANSPORT_XSECTIONS;
  G = 0;
  L = 0;

  diffusion_initialized = false;
  scattering_initialized = false;
}

//###################################################################
/**Makes a simple material with no transfer matrix just sigma_t.*/
void chi_physics::TransportCrossSections::
  MakeSimple0(int in_G, double in_sigmat)
{
  G = in_G;
  sigma_tg.resize(in_G,in_sigmat);

  transfer_matrix.push_back(chi_math::SparseMatrix(in_G,in_G));
}

//###################################################################
/**Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering.*/
void chi_physics::TransportCrossSections::
  MakeSimple1(int in_G, double in_sigmat, double c)
{
  G = in_G;

  sigma_tg.resize(in_G,in_sigmat);

  transfer_matrix.push_back(chi_math::SparseMatrix(in_G,in_G));

  if (G == 1)
    transfer_matrix[0].SetDiagonal(std::vector<double>(in_G,in_sigmat*c));
  else
    transfer_matrix[0].SetDiagonal(std::vector<double>(in_G,in_sigmat*c*2.0/4.0));

  for (int g=0; g<in_G; g++)
  {
    //Downscattering
    if (g>0)
    {
      transfer_matrix[0][g][g-1] = in_sigmat*c*2.0/4.0;
    }


    //Upscattering
    if (g>(in_G/2))
    {
      if (g<(in_G-1))
      {
        transfer_matrix[0][g][g-1] = in_sigmat*c*1.0/4.0;
        transfer_matrix[0][g][g+1] = in_sigmat*c*1.0/4.0;
      }
      else
      {
        transfer_matrix[0][g][g-1] = in_sigmat*c*2.0/4.0;
      }
    }

  }

//  chi_log.Log(LOG_0WARNING) << transfer_matrix[0].PrintS();
}