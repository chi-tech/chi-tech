#include "material_property_transportxsections.h"

#include <ChiMath/GolubFischer/GolubFischer.h>
#include <ChiMath/chi_math.h>

#include <chi_log.h>

extern ChiLog& chi_log;

#include <algorithm>

//###################################################################
/**Computes the discrete scattering tables.*/
void chi_physics::TransportCrossSections::ComputeDiscreteScattering(int in_L)
{
  if (scattering_initialized)
    return;

  chi_log.Log(LOG_0) << "Creating Discrete scattering angles.";

  //============================================= Compute scattering energy
  //                                              probabilities
  std::vector<std::vector<double>> prob_gprime_g;        //Dense version of S[0]
  std::vector<std::vector<double>> prob_gprime_g_normed; //Normalized

  prob_gprime_g.       resize(num_groups, std::vector<double>(num_groups, 0.0));
  prob_gprime_g_normed.resize(num_groups, std::vector<double>(num_groups, 0.0));
  std::vector<double>  prob_sum(num_groups, 0.0);

  std::vector<std::vector<double>> prob_g_gprime_normed; //Transpose to be computed
  prob_g_gprime_normed.resize(num_groups, std::vector<double>(num_groups, 0.0));

  //============================================= Extract the dense version
  for (int g=0; g < num_groups; g++)
  {
    int num_transfer = transfer_matrix[0].rowI_indices[g].size();
    for (int j=0; j<num_transfer; j++)
    {
      int gp = transfer_matrix[0].rowI_indices[g][j];
      prob_gprime_g[g][gp] = transfer_matrix[0].rowI_values[g][j];
    }//for j
  }//for g

  //============================================= Compute the column sum
  for (int gp=0; gp < num_groups; gp++)
  {
    for (int g=0; g < num_groups; g++)
    {
      prob_gprime_g_normed[g][gp] = prob_gprime_g[g][gp];

      if (g>0) prob_gprime_g_normed[g][gp] += prob_gprime_g_normed[g-1][gp];

      prob_sum[gp] += prob_gprime_g[g][gp];
    }//for gp
  }//for g

  //============================================= Compute normalization (CDF)
  //This step just normalizes a given row so
  //that its L1-norm is 1.0
  for (int g=0; g < num_groups; g++)
    for (int gp=0; gp < num_groups; gp++)
      prob_gprime_g_normed[g][gp] /= prob_sum[gp];


  //============================================= copy cdf
  //This will make the CDF from 0.0 to 1.0
  cdf_gprime_g.resize(num_groups, std::vector<double>(num_groups, 0.0));
  for (int gp=0; gp < num_groups; gp++)
    for (int g=0; g < num_groups; g++)
      cdf_gprime_g[g][gp] = prob_gprime_g_normed[g][gp];

   cdf_gprime_g = chi_math::Transpose(cdf_gprime_g);

  //============================================= Collect scattering moments
  //For a given scattering gprime->g we need all
  //the moments so we can construct the angles
  std::vector<std::vector<std::vector<double>>> moment_gprime_g_m;
  int max_ord = in_L;
  if (in_L > scattering_order) max_ord = scattering_order;

  moment_gprime_g_m.resize(num_groups);
  for (int gp=0; gp < num_groups; gp++)
  {
    moment_gprime_g_m[gp].resize(num_groups);
    for (int g=0; g < num_groups; g++)
    {
      moment_gprime_g_m[gp][g].resize(max_ord+1,0.0);
      for (int m=0; m<=max_ord; m++)
      {
        moment_gprime_g_m[gp][g][m] = transfer_matrix[m].ValueIJ(g,gp);
      }
    }
  }

  //============================================= Compute angles
  GolubFischer gb;
  scat_angles_gprime_g.resize(num_groups);
  for (int gp=0; gp < num_groups; gp++)
  {
    scat_angles_gprime_g[gp].resize(num_groups);

    for (int g=0; g < num_groups; g++)
    {
      if (transfer_matrix[0].ValueIJ(g,gp) < 1.0e-16) continue;

      Tvecdbl_vecdbl angles = gb.GetDiscreteScatAngles(moment_gprime_g_m[gp][g]);
      std::copy(angles.begin(),angles.end(),
                std::back_inserter(scat_angles_gprime_g[gp][g]));

      int num_angles = scat_angles_gprime_g[gp][g].size();
      double intgl = 0.0;
      for (int a=0; a<num_angles; a++)
        intgl += scat_angles_gprime_g[gp][g][a].second;


      double cumulative = 0.0;
      for (int a=0; a<num_angles; a++)
      {
        cumulative += scat_angles_gprime_g[gp][g][a].second/intgl;
        scat_angles_gprime_g[gp][g][a].second = cumulative;
      }
    }
  }

  chi_log.Log(LOG_0) << "Done creating Discrete scattering angles.";
  scattering_initialized = true;
}


//###################################################################
/**Samples the g to gprime scattering probability to determine exit
 * energy.*/
int chi_physics::TransportCrossSections::Sample_gprime(int gp, double rn)
{
  int gto = 0;

  gto = std::lower_bound(cdf_gprime_g[gp].begin(),
                         cdf_gprime_g[gp].end(),
                         rn) - cdf_gprime_g[gp].begin();

  return gto;
}

//###################################################################
/**Sample scattering angles.*/
double chi_physics::TransportCrossSections::
  SampleMu_gprime_g(int gp, int g, double rn, bool isotropic)
{
  double mu = 0.0;

  struct
  {
    bool operator()(const std::pair<double,double>& left, double val)
    {return left.second <= val;}
  }compare;

  if (isotropic or scat_angles_gprime_g[gp][g].empty())
    mu = 2.0*rn-1.0;
  else
  {
    int angle_num = std::lower_bound(scat_angles_gprime_g[gp][g].begin(),
                                     scat_angles_gprime_g[gp][g].end(),rn,
                                     compare) -
                                     scat_angles_gprime_g[gp][g].begin();
    mu = scat_angles_gprime_g[gp][g][angle_num].first;
  }

  if (std::isnan(mu))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Mu corruption in GolubFischer sample mu.\n"
      << " g      =" << g << "\n"
      << " gprime =" << gp << "\n"
      << " rn     =" << rn;
  }

  if (mu > 1.0)
    mu = 1.0;
  if (mu < -1.0)
    mu = -1.0;


  return mu;
}

