#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"


//###################################################################
/**Partial Jacobi energy collapse.*/
void chi_physics::TransportCrossSections::
  EnergyCollapse(std::vector<double>& ref_xi,
                 double& D, double& ref_sigma_a,
                 int collapse_type)
{
  //============================================= Make a Dense matrix from
  //                                              sparse transfer matrix
  std::vector<std::vector<double>> S;
  S.resize(num_groups, std::vector<double>(num_groups, 0.0));
  for (int g=0; g < num_groups; g++)
  {
    S[g][g] = 0.0;
    size_t num_transfer = transfer_matrices[0].rowI_indices[g].size();
    for (size_t j=0; j<num_transfer; j++)
    {
      size_t gprime   = transfer_matrices[0].rowI_indices[g][j];
      S[g][gprime] = transfer_matrices[0].rowI_values[g][j];
    }//for j
  }//for g

  //============================================= Compiling the A and B matrices
  //                                              for different methods
  MatDbl A(num_groups, VecDbl(num_groups, 0.0));
  MatDbl B(num_groups, VecDbl(num_groups, 0.0));
  for (int g=0; g < num_groups; g++)
  {
    if      (collapse_type == E_COLLAPSE_JACOBI)
    {
//      A[g][g] = sigma_t[g] - S[g][g];
//      for (int gp=0; gp<g; gp++)
//        B[g][gp] = S[g][gp];
//
//      for (int gp=g+1; gp < num_groups; gp++)
//        B[g][gp] = S[g][gp];

      A[g][g] = sigma_t[g];
      B = S;
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_JACOBI)
    {
      A[g][g] = sigma_t[g];
      for (int gp=0; gp < num_groups; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_GAUSS)
    {
      A[g][g] = sigma_t[g] - S[g][g];
      for (int gp=0; gp<g; gp++)
        A[g][gp] = -S[g][gp];

      for (int gp=g+1; gp < num_groups; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_GAUSS)
    {
      A[g][g] = sigma_t[g];
      for (int gp=0; gp<g; gp++)
        A[g][gp] = -S[g][gp];

      for (int gp=g; gp < num_groups; gp++)
        B[g][gp] = S[g][gp];
    }
  }//for g

  //============================================= Correction for zero xs groups
  //Some cross-sections developed from monte-carlo
  //methods can result in some of the groups
  //having zero cross-sections. In that case
  //it will screw up the power iteration
  //initial guess of 1.0. Here we reset them
  for (int g=0; g < num_groups; g++)
    if (sigma_t[g] < 1.0e-16)
      A[g][g] = 1.0;

  MatDbl Ainv = chi_math::Inverse(A);
  MatDbl C    = chi_math::MatMul(Ainv,B);
  VecDbl E(num_groups, 1.0);

  //============================================= Perform power iteration
  double rho = chi_math::PowerIteration(C, E, 1000, 1.0e-12);

  ref_xi.resize(num_groups, 0.0);
  double sum = 0.0;
  for (int g=0; g < num_groups; g++)
    sum += std::fabs(E[g]);

  for (int g=0; g < num_groups; g++)
    ref_xi[g] = std::fabs(E[g])/sum;


  //======================================== Compute two-grid diffusion quantities
  D = 0.0;
  ref_sigma_a = 0.0;
  for (int g=0; g < num_groups; ++g)
  {
    D += diffusion_coeff[g] * ref_xi[g];

    ref_sigma_a += sigma_t[g] * ref_xi[g];

    for (int gp=0; gp < num_groups; ++gp)
      ref_sigma_a -= S[g][gp] * ref_xi[gp];
  }

  //======================================== Verbose output the spectrum
  chi::log.Log0Verbose1() << "Fundamental eigen-value: " << rho;
  std::stringstream outstr;
  for (auto& xi : ref_xi)
    outstr << xi << '\n';
  chi::log.Log0Verbose1() << outstr.str();
}
