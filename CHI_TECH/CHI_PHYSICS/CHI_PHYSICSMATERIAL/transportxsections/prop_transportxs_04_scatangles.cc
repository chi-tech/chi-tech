#include "../property10_transportxsections.h"

#include <CHI_MATH/GolubFischer/GolubFischer.h>
#include <CHI_MATH/chi_math.h>

#include <chi_log.h>

extern CHI_LOG chi_log;

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

  prob_gprime_g.       resize(G,std::vector<double>(G,0.0));
  prob_gprime_g_normed.resize(G,std::vector<double>(G,0.0));
  std::vector<double>  prob_sum(G,0.0);

  std::vector<std::vector<double>> prob_g_gprime_normed; //Transpose to be computed
  prob_g_gprime_normed.resize(G,std::vector<double>(G,0.0));

  //============================================= Extract the dense version
  for (int g=0; g<G; g++)
  {
    int num_transfer = transfer_matrix[0].inds_rowI[g].size();
    for (int j=0; j<num_transfer; j++)
    {
      int gp = transfer_matrix[0].inds_rowI[g][j];
      prob_gprime_g[g][gp] = transfer_matrix[0].rowI_colJ[g][j];
    }//for j
  }//for g

  //============================================= Compute the column sum
  for (int gp=0; gp<G; gp++)
  {
    for (int g=0; g<G; g++)
    {
      prob_gprime_g_normed[g][gp] = prob_gprime_g[g][gp];

      if (g>0) prob_gprime_g_normed[g][gp] += prob_gprime_g_normed[g-1][gp];

      prob_sum[gp] += prob_gprime_g[g][gp];
    }//for gp
  }//for g

  //============================================= Compute normalization (CDF)
  //This step just normalizes a given row so
  //that its L1-norm is 1.0
  for (int g=0; g<G; g++)
    for (int gp=0; gp<G; gp++)
      prob_gprime_g_normed[g][gp] /= prob_sum[gp];


  //============================================= copy cdf
  //This will make the CDF from 0.0 to 1.0
  cdf_gprime_g.resize(G,std::vector<double>(G,0.0));
  double intgl = 0.0;
  for (int gp=0; gp<G; gp++)
    for (int g=0; g<G; g++)
    {
      cdf_gprime_g[g][gp] = prob_gprime_g_normed[g][gp];
//      if (g<62 and gp<62)
//      {
//        chi_log.Log(LOG_0)
//        << "gp=" << gp << " g=" << g << " cp=" << cdf_gprime_g[g][gp];
//      }
    }



  //============================================= Collect scattering moments
  //For a given scattering gprime->g we need all
  //the moments so we can construct the angles
  std::vector<std::vector<std::vector<double>>> moment_gprime_g_m;
  int max_ord = in_L;
  if (in_L > L) max_ord = L;

  moment_gprime_g_m.resize(G);
  for (int gp=0; gp<G; gp++)
  {
    moment_gprime_g_m[gp].resize(G);
    for (int g=0; g<G; g++)
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
  scat_angles_gprime_g.resize(G);
  for (int gp=0; gp<G; gp++)
  {
    scat_angles_gprime_g[gp].resize(G);
//    if (sigma_tg[gp]<1.0e-16) continue;

    for (int g=0; g<G; g++)
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
  for (int g=0; g<G; g++)
  {
    if ((rn < cdf_gprime_g[g][gp]) and (gp == 0))
    {gto = g; break;}

    if ((rn < cdf_gprime_g[g][gp]) and (rn >= cdf_gprime_g[g-1][gp]))
    {gto = g; break;}
  }

  return gto;
}

//###################################################################
/**Sample scattering angles.*/
double chi_physics::TransportCrossSections::
  SampleMu_gprime_g(int gp, int g, double rn, bool isotropic)
{
  double mu = 0.0;

  if (isotropic or scat_angles_gprime_g[gp][g].size()==0)
    mu = 2.0*rn-1.0;
  else
  {
    int num_angles = scat_angles_gprime_g[gp][g].size();
    for (int a=0; a<num_angles; a++)
    {
      double prob = scat_angles_gprime_g[gp][g][a].second;
      if ((rn < prob) and (a == 0))
      {
        mu =  scat_angles_gprime_g[gp][g][a].first;
        break;
      }
      if ((rn < prob) and (rn >= scat_angles_gprime_g[gp][g][a-1].second))
      {
        mu =  scat_angles_gprime_g[gp][g][a].first;
        break;
      }
    }
  }

  if (isnan(mu))
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

