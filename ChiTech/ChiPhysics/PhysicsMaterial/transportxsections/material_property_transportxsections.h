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

  //Monte-Carlo quantities
public:
  bool scattering_initialized = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

public:
  //00
  TransportCrossSections();

private:
  void Reset();
  void ComputeAbsorption();

  public:
  void MakeSimple0(int in_G, double in_sigmat);
  void MakeSimple1(int in_G, double in_sigmat, double c);
  void MakeCombined(std::vector<std::pair<int,double>>& combinations);

  //01
  void MakeFromPDTxsFile(const std::string &file_name,const std::string& MT_TRANSFER);
  void MakeFromCHIxsFile(const std::string &file_name);
  void FinalizeCrossSections();

  //02
  void ComputeDiffusionParameters();

  //05
  void PushLuaTable(lua_State* L) override;

  //06
  void ExportToChiFormat(const std::string& file_name);


};

}//namespace chi_physics


#endif