#ifndef CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H
#define CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H

#include "ChiPhysics/PhysicsMaterial/material_property_base.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;

/**\defgroup LuaTransportXSs Transport Cross Sections
 * \ingroup LuaPhysics*/

namespace chi_physics
{

//###################################################################
/** Class for handling Transport-Theory related cross sections.*/
class TransportCrossSections : public chi_physics::MaterialProperty
{
protected:
  using TransferMatrix = chi_math::SparseMatrix;
  using EmissionSpectra = std::vector<std::vector<double>>;

public:
  /// A struct containing data for a delayed neutron precursor.
  struct Precursor
  {
    double decay_constant = 0.0;
    double fractional_yield = 0.0;
    std::vector<double> emission_spectrum;
  };

public:
  unsigned int num_groups_ = 0;       ///< Total number of groups
  unsigned int scattering_order_ = 0; ///< Legendre scattering order
  unsigned int num_precursors_ = 0;   ///< Number of precursors

  bool is_fissionable_ = false;
  bool is_fission_scaled_ = false;

  /// Energy bin boundaries in MeV
  std::vector<std::vector<double>> e_bounds_;

  std::vector<double> sigma_t_;  ///< Total cross section
  std::vector<double> sigma_a_;  ///< Absorption cross section
  std::vector<double> sigma_f_;  ///< Fission cross section

  std::vector<double> nu_sigma_f_;
  std::vector<double> nu_prompt_sigma_f_;
  std::vector<double> nu_delayed_sigma_f_;

  std::vector<double> inv_velocity_;

  std::vector<TransferMatrix> transfer_matrices_;
  std::vector<std::vector<double>> production_matrix_;

  std::vector<Precursor> precursors_;

  //Diffusion quantities
public:
  bool diffusion_initialized_ = false;
public:
  std::vector<double> diffusion_coeff_; ///< Transport corrected diffusion coeff
  std::vector<double> sigma_removal_;   ///< Removal cross section
  std::vector<double> sigma_s_gtog_;    ///< Within-group scattering xs

  //Monte-Carlo quantities
public:
  bool scattering_initialized_ = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g_;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g_;

public:
  //00
  TransportCrossSections();

private:
  void Reset();

public:
  void MakeSimple0(int n_grps, double sigma);
  void MakeSimple1(int n_grps, double sigma, double c);
  void MakeCombined(std::vector<std::pair<int,double>>& combinations);

private:
  void ComputeAbsorption();

public:
  void ScaleFissionData(double k);

public:
  //01
  void MakeFromChiXSFile(const std::string &file_name);

  //02
  void ComputeDiffusionParameters();

  //05
  void PushLuaTable(lua_State* L) override;

  //06
  void ExportToChiXSFile(const std::string& file_name);


};

}//namespace chi_physics


#endif