#include "../property10_transportxsections.h"

#include <CHI_MATH/../../CHI_RESOURCES/Dependencies/Eigen/Dense>

#include <chi_log.h>

extern CHI_LOG chi_log;

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
    int num_transfer = transfer_matrix[0].inds_rowI[g].size();
    for (int j=0; j<num_transfer; j++)
    {
      int gprime   = transfer_matrix[0].inds_rowI[g][j];
      S[g][gprime] = transfer_matrix[0].rowI_colJ[g][j];
    }//for j
  }//for g

  //============================================= Compiling the A and B matrices
  //                                              for different methods
  Eigen::MatrixXd A(G,G); A.setZero();
  Eigen::MatrixXd B(G,G); B.setZero();
  for (int g=0; g<G; g++)
  {
    if      (collapse_type == E_COLLAPSE_JACOBI)
    {
      A(g,g) = sigma_tg[g] - S[g][g];
      for (int gp=0; gp<g; gp++)
        B(g,gp) = S[g][gp];

      for (int gp=g+1; gp<G; gp++)
        B(g,gp) = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_JACOBI)
    {
      A(g,g) = sigma_tg[g];
      for (int gp=0; gp<G; gp++)
        B(g,gp) = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_GAUSS)
    {
      A(g,g) = sigma_tg[g] - S[g][g];
      for (int gp=0; gp<g; gp++)
        A(g,gp) = -S[g][gp];

      for (int gp=g+1; gp<G; gp++)
        B(g,gp) = S[g][gp];
    }
    else if (collapse_type == E_COLLAPSE_PARTIAL_GAUSS)
    {
      A(g,g) = sigma_tg[g];
      for (int gp=0; gp<g; gp++)
        A(g,gp) = -S[g][gp];

      for (int gp=g; gp<G; gp++)
        B(g,gp) = S[g][gp];
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
      A(g,g) = 1.0;

  Eigen::MatrixXd Ainv = A.inverse();
  Eigen::MatrixXd C = Ainv*B;
  Eigen::VectorXd E(G); E.setConstant(1.0);
  Eigen::VectorXd E_new(G);

  //============================================= Perform power iteration
  E_new = C*E;
  double rho = E_new.dot(E);
  double rho_old;
  E = E_new/E_new.norm();

  if (rho<0.0)
    E = E*-1.0;

  for (int k=0; k<1000; k++)
  {
    rho_old = rho;

    E_new = C*E;
    rho = E_new.dot(E);
    E = E_new/E_new.norm();

    if (rho<0.0)
      E = E*-1.0;

    if (std::fabs(rho-rho_old) < 1.0e-12)
      break;
  }
  E = E/rho;

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
}