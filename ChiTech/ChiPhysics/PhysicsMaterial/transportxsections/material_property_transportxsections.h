#ifndef CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H
#define CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H

#include "ChiPhysics/PhysicsMaterial/material_property_base.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#define E_COLLAPSE_PARTIAL_JACOBI 1
#define E_COLLAPSE_JACOBI         2
#define E_COLLAPSE_PARTIAL_GAUSS  3
#define E_COLLAPSE_GAUSS          4

typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;

namespace chi_physics
{

//###################################################################
/** Basic thermal conductivity material property.*/
class TransportCrossSections : public chi_physics::MaterialProperty
{
public:
  size_t num_groups=0;                     ///< Total number of Groups
  size_t scattering_order=0;               ///< Legendre scattering order
  size_t num_precursors=0;                 ///< Number of precursors
  bool is_fissile = false;                 ///< Fissile or not

  std::vector<double> sigma_tg;           ///< Total cross-section
  std::vector<double> sigma_fg;           ///< Sigmaf cross-section
  std::vector<double> sigma_captg;        ///< Capture cross-section
  std::vector<double> chi_g;              ///< Fission spectrum
  std::vector<double> nu_sigma_fg;        ///< Nubar-Sigmaf cross-section
  std::vector<double> nu_p_sigma_fg;      ///< Prompt-Nubar-Sigmaf cross-section
  std::vector<double> nu_d_sigma_fg;      ///< Delayed-Nubar-Sigmaf cross-section
  std::vector<double> ddt_coeff;          ///< Time derivative coefficient
  std::vector<std::vector<double>> chi_d; ///< Delayed neutron fission spectrum
  std::vector<double> lambda;             ///< Delayed neutron decay constants
  std::vector<double> gamma;              ///< Delayed neutron yields

  std::vector<chi_math::SparseMatrix> transfer_matrix;

  //Diffusion quantities
public:
  bool diffusion_initialized = false;
public:
  std::vector<double> diffg;               ///< Transport corrected Diffusion coeff
  std::vector<double> sigma_rg;            ///< Removal cross-section
  std::vector<double> sigma_ag;            ///< Pure absorption
  std::vector<double> sigma_s_gtog;        ///< Within-group scattering xs

  //Two-grid acceleration quantities
  std::vector<double> xi_Jfull_g;          ///< Infinite medium spectrum Jfull
  std::vector<double> xi_Jpart_g;          ///< Infinite medium spectrum Jpartial

  double D_jfull = 0.0;                    ///< Collapsed Diffusion coefficient Jfull
  double D_jpart = 0.0;                    ///< Collapsed Diffusion coefficient Jpart

  double sigma_a_jfull = 0.0;              ///< Collapsed absorption Jfull
  double sigma_a_jpart = 0.0;              ///< Collapsed absorption Jpart

  //Monte-Carlo quantities
public:
  bool scattering_initialized = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

private:
  void Reset()
  {
    num_groups= scattering_order= num_precursors=0;

    sigma_tg.clear();
    sigma_fg = sigma_captg = chi_g = nu_sigma_fg = ddt_coeff = sigma_tg;
    nu_p_sigma_fg = nu_d_sigma_fg = nu_sigma_fg;
    transfer_matrix.clear();
    lambda.clear();
    gamma.clear();
    chi_d.clear();

    diffusion_initialized = false;
    scattering_initialized = false;

    diffg = sigma_rg = sigma_ag = sigma_s_gtog = xi_Jfull_g = xi_Jpart_g = sigma_tg;

    D_jfull = D_jpart = sigma_a_jfull = sigma_a_jpart = 0.0;

    cdf_gprime_g.clear();
    scat_angles_gprime_g.clear();
  }
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

}//namespace chi_physics


#endif