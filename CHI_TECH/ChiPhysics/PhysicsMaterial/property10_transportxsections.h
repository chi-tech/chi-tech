#ifndef _chi_physics_transportxsections_h
#define _chi_physics_transportxsections_h

#include "chi_physicsmaterial.h"
#include <ChiMath/SparseMatrix/chi_math_sparse_matrix.h>

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
  int G=0;
  int L=0;

  std::vector<double> sigma_tg;     ///< MT 1    Total cross-section
  std::vector<double> sigma_fg;     ///< MT 18   Sigmaf cross-section
  std::vector<double> sigma_captg;  ///< MT 27   Capture cross-section
  std::vector<double> chi_g;        ///< MT 2018 Fission spectrum
  std::vector<double> nu_sigma_fg;  ///< MT 2452 Nubar-Sigmaf cross-section
  std::vector<double> ddt_coeff;    ///< Time derivative coefficient

  /**The MT number for this transfer varies:
   * MT 2500 is total,
   * MT 2501 is scattering only
   * MT 2519 is scattering and fission
   * MT 2502 is elastic scattering only
   * MT 2504 is inelastic scattering only*/
  std::vector<chi_math::SparseMatrix> transfer_matrix;

public:
  bool diffusion_initialized = false;
  bool scattering_initialized = false;
public:
  std::vector<double> diffg;
  std::vector<double> sigma_rg;
  std::vector<double> sigma_ag;
  std::vector<double> sigma_s_gtog;


  std::vector<double> xi_Jfull_g;
  std::vector<double> xi_Jpart_g;

  double D_jfull = 0.0;
  double D_jpart = 0.0;

  double sigma_a_jfull = 0.0;
  double sigma_a_jpart = 0.0;

private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

public:
  //00
  TransportCrossSections();


  void MakeSimple0(int in_G, double in_sigmat);
  void MakeSimple1(int in_G, double in_sigmat, double c);
  void MakeCombined(std::vector<std::pair<int,double>>& combinations);

  //01
  void MakeFromPDTxsFile(const std::string &file_name,std::string MT_TRANSFER);
  void MakeFromCHIxsFile(const std::string &file_name);

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

  //05
  void PushLuaTable(lua_State* L) override;


};


#endif