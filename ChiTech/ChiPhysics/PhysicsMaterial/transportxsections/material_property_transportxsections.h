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
/** Class for handling Transport-Theory related cross sections.*/
class TransportCrossSections : public chi_physics::MaterialProperty
{
public:
  size_t num_groups=0;                          ///< Total number of Groups
  size_t scattering_order=0;                    ///< Legendre scattering order
  size_t num_precursors=0;                      ///< Number of precursors
  bool is_fissile = false;                      ///< Fissile or not

  typedef double GrpVal;                        ///< Denoting value per group
  typedef double PrecursorVal;                  ///< Denoting value per precursor

  std::vector<GrpVal> sigma_t;                  ///< Total cross section
  std::vector<GrpVal> sigma_f;                  ///< Sigmaf cross section
  std::vector<GrpVal> sigma_a;                  ///< Pure absorption
  std::vector<GrpVal> chi;                      ///< Fission spectrum
  std::vector<GrpVal> chi_prompt;               ///< Prompt fission spectrum
  std::vector<GrpVal> nu;                       ///< Nubar
  std::vector<GrpVal> nu_prompt;                ///< Nubar-prompt
  std::vector<GrpVal> nu_delayed;               ///< Nubar-delayed
  std::vector<GrpVal> nu_sigma_f;               ///< Nubar-Sigmaf cross section
  std::vector<GrpVal> nu_prompt_sigma_f;        ///< Prompt-Nubar-Sigmaf cross section
  std::vector<GrpVal> nu_delayed_sigma_f;       ///< Delayed-Nubar-Sigmaf cross section
  std::vector<GrpVal> inv_velocity;             ///< Groupwise inverse velocities

  std::vector<chi_math::SparseMatrix> transfer_matrices;

  std::vector<PrecursorVal> precursor_lambda;         ///< Delayed neutron decay constants
  std::vector<PrecursorVal> precursor_yield;          ///< Delayed neutron yields
  std::vector<std::vector<PrecursorVal>> chi_delayed; ///< Delayed neutron fission spectrum

  //Diffusion quantities
public:
  bool diffusion_initialized = false;
public:
  std::vector<GrpVal> diffusion_coeff; ///< Transport corrected Diffusion coeff
  std::vector<GrpVal> sigma_removal;   ///< Removal cross section
  std::vector<GrpVal> sigma_s_gtog;    ///< Within-group scattering xs

  //Two-grid acceleration quantities
  std::vector<GrpVal> xi_Jfull;        ///< Infinite medium spectrum full Jacobi-Splitting
  std::vector<GrpVal> xi_Jpart;        ///< Infinite medium spectrum partial Jacobi-Splitting

  double D_jfull = 0.0;                ///< Collapsed Diffusion coeff full Jacobi-Splitting
  double D_jpart = 0.0;                ///< Collapsed Diffusion coeff partial Jacobi-Splitting

  double sigma_a_jfull = 0.0;          ///< Collapsed absorption full Jacobi-Splitting
  double sigma_a_jpart = 0.0;          ///< Collapsed absorption partial Jacobi-Splitting

  //Monte-Carlo quantities
public:
  bool scattering_initialized = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

private:
  void Reset()
  {
    num_groups = 0;
    scattering_order = 0;
    num_precursors = 0;
    is_fissile = false;

    sigma_t.clear();
    sigma_f.clear();
    sigma_a.clear();
    chi.clear();
    chi_prompt.clear();
    nu.clear();
    nu_prompt.clear();
    nu_delayed.clear();
    nu_sigma_f.clear();
    nu_prompt_sigma_f.clear();
    nu_delayed_sigma_f.clear();
    inv_velocity.clear();
    precursor_lambda.clear();
    precursor_yield.clear();
    chi_delayed.clear();

    transfer_matrices.clear();

    //Diffusion quantities
    diffusion_initialized = false;
    diffusion_coeff.clear();
    sigma_removal.clear();
    sigma_s_gtog.clear();

    //Two-grid acceleration quantities
    xi_Jfull.clear();
    xi_Jpart.clear();

    D_jfull = 0.0;
    D_jpart = 0.0;

    sigma_a_jfull = 0.0;
    sigma_a_jpart = 0.0;

    //Monte-Carlo quantities
    scattering_initialized = false;
    cdf_gprime_g.clear();
    scat_angles_gprime_g.clear();
  }

  std::vector<GrpVal> ComputeAbsorptionXSFromTransfer();

public:
  //00
  TransportCrossSections();

  void MakeSimple0(int in_G, double in_sigmat);
  void MakeSimple1(int in_G, double in_sigmat, double c);
  void MakeCombined(std::vector<std::pair<int,double>>& combinations);

  //01
  void MakeFromPDTxsFile(const std::string &file_name,const std::string& MT_TRANSFER);
  void MakeFromCHIxsFile(const std::string &file_name);
  void FinalizeCrossSections();

  //02
  void ComputeDiffusionParameters();

  //03
  void EnergyCollapse(std::vector<double>& ref_xi,
                      double& D, double& sigma_a,
                      int collapse_type = E_COLLAPSE_JACOBI);

  //05
  void PushLuaTable(lua_State* L) override;

  //06
  void ExportToChiFormat(const std::string& file_name);


};

}//namespace chi_physics


#endif