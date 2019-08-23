#ifndef _chi_physics_transportxsections_h
#define _chi_physics_transportxsections_h

#include "chi_physicsmaterial.h"
#include <CHI_MATH/SparseMatrix/chi_math_sparse_matrix.h>

#define E_COLLAPSE_PARTIAL_JACOBI 1
#define E_COLLAPSE_JACOBI         2
#define E_COLLAPSE_PARTIAL_GAUSS  3
#define E_COLLAPSE_GAUSS          4

typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;

//###################################################################
/** Basic thermal conductivity material property.*/
class chi_physics::TransportCrossSections : public chi_physics::MaterialProperty
{
public:
  int G;
  int L;

  std::vector<double> sigma_tg;
  std::vector<chi_math::SparseMatrix> transfer_matrix;

public:
  bool diffusion_initialized;
  bool scattering_initialized;
public:
  std::vector<double> diffg;
  std::vector<double> sigma_rg;
  std::vector<double> sigma_ag;
  std::vector<double> sigma_s_gtog;

  std::vector<double> xi_Jfull_g;
  std::vector<double> xi_Jpart_g;

  double D_jfull;
  double D_jpart;

  double sigma_a_jfull;
  double sigma_a_jpart;

private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

public:
  //00
  TransportCrossSections();


  void MakeSimple0(int in_G, double in_sigmat);
  void MakeSimple1(int in_G, double in_sigmat, double c);

  //01
  void MakeFromPDTxsFile(const std::string &file_name,std::string MT_TRANSFER);

  //02
  void ComputeDiffusionParameters();

  //03
  void EnergyCollapse(std::vector<double>& ref_xi,
                      double& D, double& sigma_a,
                      int collapse_type = E_COLLAPSE_JACOBI);

  //04
  void ComputeDiscreteScattering(int in_L);
  int  Sample_gprime(int g,double rn);
  double SampleMu_gprime_g(int g, int gprime, double rn, bool isotropic = false);
};


#endif