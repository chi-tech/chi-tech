#include "../property10_transportxsections.h"

#include <chi_log.h>

extern ChiLog chi_log;


//###################################################################
/**Partial Jacobi energy collapse.*/
void chi_physics::TransportCrossSections::
  EnergyCollapse(std::vector<double>& ref_xi,
                 double& D, double& sigma_a,
                 int collapse_type)
{
  //============================================= Make a Dense matrix from
  //                                              sparse transfer matrix
  std::vector<std::vector<double>> S;
  S.resize(G,std::vector<double>(G,0.0));
  for (int g=0; g<G; g++)
  {
    S[g][g] = 1.0;
    int num_transfer = transfer_matrix[0].rowI_indices[g].size();
    for (int j=0; j<num_transfer; j++)
    {
      int gprime   = transfer_matrix[0].rowI_indices[g][j];
      S[g][gprime] = transfer_matrix[0].rowI_values[g][j];
    }//for j
  }//for g

  //============================================= Compiling the A and B matrices
  //                                              for different methods
  MatDbl A(G, VecDbl(G,0.0));
  MatDbl B(G, VecDbl(G,0.0));
  for (int g=0; g<G; g++)
  {
    if      (collapse_type == E_COLLAPSE_JACOBI)
    {
      A[g][g] = sigma_tg[g] - S[g][g];
      for (int gp=0; gp<g; gp++)
        B[g][gp] = S[g][gp];

      for (int gp=g+1; gp<G; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_JACOBI)
    {
      A[g][g] = sigma_tg[g];
      for (int gp=0; gp<G; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_GAUSS)
    {
      A[g][g] = sigma_tg[g] - S[g][g];
      for (int gp=0; gp<g; gp++)
        A[g][gp] = -S[g][gp];

      for (int gp=g+1; gp<G; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_GAUSS)
    {
      A[g][g] = sigma_tg[g];
      for (int gp=0; gp<g; gp++)
        A[g][gp] = -S[g][gp];

      for (int gp=g; gp<G; gp++)
        B[g][gp] = S[g][gp];
    }
  }//for g

  //============================================= Correction for zero xs groups
  //Some cross-sections developed from monte-carlo
  //methods can result in some of the groups
  //having zero cross-sections. In that case
  //it will screw up the power iteration
  //initial guess of 1.0. Here we reset them
  for (int g=0; g<G; g++)
    if (sigma_tg[g] < 1.0e-16)
      A[g][g] = 1.0;

  MatDbl Ainv = chi_math::Inverse(A);
  MatDbl C    = chi_math::MatMul(Ainv,B);
  VecDbl E(G,1.0);

  //============================================= Perform power iteration
  double rho = chi_math::PowerIteration(C, E, 1000, 1.0e-12);

  ref_xi.resize(G);
  double sum = 0.0;
  for (int g=0; g<G; g++)
    sum += E[g];

  for (int g=0; g<G; g++)
    ref_xi[g] = E[g]/sum;


  //======================================== Compute two-grid diffusion quantities
  D = 0.0;
  sigma_a = 0.0;
  for (int g=0; g<G; g++)
  {
    D += diffg[g]*ref_xi[g];

    sigma_a += sigma_tg[g]*ref_xi[g];

    for (int gp=0; gp<G; gp++)
      sigma_a -= S[g][gp]*ref_xi[gp];
  }

  //======================================== Verbose output the spectrum
  chi_log.Log(LOG_0VERBOSE_1) << "Fundamental eigen-value: " << rho;
  std::stringstream outstr;
  for (auto& xi : ref_xi)
    outstr << xi << '\n';
  chi_log.Log(LOG_0VERBOSE_1) << outstr.str();

}
