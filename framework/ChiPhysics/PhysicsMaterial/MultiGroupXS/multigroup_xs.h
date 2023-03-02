#ifndef CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H
#define CHI_PHYSICS_TRANSPORT_CROSS_SECTIONS_H

#include "ChiPhysics/PhysicsMaterial/material_property_base.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"


/**\defgroup LuaTransportXSs Transport Cross Sections
 * \ingroup LuaPhysics*/

namespace chi_physics
{

//######################################################################
class MultiGroupXSBase : public MaterialProperty
{
public:
  /**
   * A struct containing data for a delayed neutron precursor.
   */
  struct Precursor
  {
    double decay_constant = 0.0;
    double fractional_yield = 0.0;
    std::vector<double> emission_spectrum;
  };

  MultiGroupXSBase() : MaterialProperty(PropertyType::TRANSPORT_XSECTIONS) {}

  virtual const unsigned int NumGroups() const = 0;
  virtual const unsigned int ScatteringOrder() const = 0;
  virtual const unsigned int NumPrecursors() const = 0;

  virtual const bool IsFissionable() const = 0;
  virtual const bool IsFissionScaled() const = 0;
  virtual const bool DiffusionInitialized() const  = 0;
  virtual const bool ScatteringInitialized() const = 0;

  virtual const std::vector<double>& SigmaTotal() const = 0;
  virtual const std::vector<double>& SigmaAbsorption() const = 0;
  virtual const std::vector<double>& SigmaFission() const = 0;

  virtual const std::vector<double>& NuSigmaF() const = 0;
  virtual const std::vector<double>& NuPromptSigmaF() const = 0;
  virtual const std::vector<double>& NuDelayedSigmaF() const = 0;

  virtual const std::vector<double>& InverseVelocity() const = 0;

  virtual const std::vector<chi_math::SparseMatrix>&
  TransferMatrices() const = 0;

  virtual const chi_math::SparseMatrix&
  TransferMatrix(unsigned int ell) const = 0;

  virtual const std::vector<std::vector<double>> ProductionMatrix() const = 0;

  virtual const std::vector<Precursor>& Precursors() const = 0;

  virtual const std::vector<double>& DiffusionCoefficient() const = 0;
  virtual const std::vector<double>& SigmaRemoval() const = 0;
  virtual const std::vector<double>& SigmaSGtoG() const = 0;
};


//######################################################################
/** Class for handling Transport-Theory related cross sections.*/
class MultiGroupXS : public MultiGroupXSBase
{
protected:
  typedef std::vector<std::pair<double,double>> AnglePairs;

protected:
  unsigned int num_groups_ = 0;       ///< Total number of groups
  unsigned int scattering_order_ = 0; ///< Legendre scattering order
  unsigned int num_precursors_ = 0;   ///< Number of precursors

  bool is_fissionable_ = false;
  bool is_fission_scaled_ = false;

  std::vector<std::vector<double>> e_bounds_; ///< Energy bin boundaries in MeV

  std::vector<double> sigma_t_;  ///< Total cross section
  std::vector<double> sigma_a_;  ///< Absorption cross section
  std::vector<double> sigma_f_;  ///< Fission cross section

  std::vector<double> nu_sigma_f_;
  std::vector<double> nu_prompt_sigma_f_;
  std::vector<double> nu_delayed_sigma_f_;

  std::vector<double> inv_velocity_;

  std::vector<chi_math::SparseMatrix> transfer_matrices_;
  std::vector<std::vector<double>> production_matrix_;

  std::vector<Precursor> precursors_;

  //Diffusion quantities
  bool diffusion_initialized_ = false;
  std::vector<double> diffusion_coeff_; ///< Transport corrected diffusion coeff
  std::vector<double> sigma_removal_;   ///< Removal cross section
  std::vector<double> sigma_s_gtog_;    ///< Within-group scattering xs

  //Monte-Carlo quantities
protected:
  bool scattering_initialized_ = false;
private:
  std::vector<std::vector<double>>         cdf_gprime_g_;
  std::vector<std::vector<AnglePairs>> scat_angles_gprime_g_;

public:
  //00
  MultiGroupXS() :
      MultiGroupXSBase(),
      num_groups_(0), scattering_order_(0), num_precursors_(0),
      diffusion_initialized_(false), scattering_initialized_(false)
  {}

private:
  void Reset();

public:
  void MakeSimple0(unsigned int num_groups, double sigma_t);
  void MakeSimple1(unsigned int num_groups, double sigma_t, double c);
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

  //Accessors
  const unsigned int NumGroups() const override { return num_groups_; }

  const unsigned int ScatteringOrder() const override
  { return scattering_order_; }

  const unsigned int NumPrecursors() const override { return num_precursors_; }

  const bool IsFissionable() const override { return is_fissionable_; }
  const bool IsFissionScaled() const override { return is_fission_scaled_; }

  const bool DiffusionInitialized() const override
  { return diffusion_initialized_; }

  const bool ScatteringInitialized() const override
  { return scattering_initialized_; }

  const std::vector<double>& SigmaTotal() const override { return sigma_t_; }
  const std::vector<double>& SigmaAbsorption() const override { return sigma_a_; }
  const std::vector<double>& SigmaFission() const override { return sigma_f_; }

  const std::vector<double>& NuSigmaF() const override { return nu_sigma_f_; }

  const std::vector<double>& NuPromptSigmaF() const override
  { return nu_prompt_sigma_f_; }

  const std::vector<double>& NuDelayedSigmaF() const override
  { return nu_delayed_sigma_f_; }

  const std::vector<double>& InverseVelocity() const override
  { return inv_velocity_; }

  const std::vector<chi_math::SparseMatrix>& TransferMatrices() const override
  { return transfer_matrices_; }

  const chi_math::SparseMatrix& TransferMatrix(unsigned int ell) const override
  { return transfer_matrices_.at(ell); }

  const std::vector<std::vector<double>> ProductionMatrix() const override
  { return production_matrix_; }

  const std::vector<Precursor>& Precursors() const override
  { return precursors_; }

  const std::vector<double>& DiffusionCoefficient() const override
  { return diffusion_coeff_; }

  const std::vector<double>& SigmaRemoval() const override
  { return sigma_removal_; }

  const std::vector<double>& SigmaSGtoG() const override
  { return sigma_s_gtog_; }
};


//######################################################################
class AdjointMultiGroupXS : public MultiGroupXSBase
{
private:
  const MultiGroupXS& xs_;
  std::vector<chi_math::SparseMatrix> transposed_transfer_matrices_;
  std::vector<std::vector<double>> transposed_production_matrices_;

public:
  explicit AdjointMultiGroupXS(const MultiGroupXS& xs);

  //Accessors
  const unsigned int NumGroups() const override { return xs_.NumGroups(); }

  const unsigned int ScatteringOrder() const override
  { return xs_.ScatteringOrder(); }

  const unsigned int NumPrecursors() const override
  { return xs_.NumPrecursors(); }

  const bool IsFissionable() const override { return xs_.IsFissionable(); }
  const bool IsFissionScaled() const override { return xs_.IsFissionScaled(); }

  const bool DiffusionInitialized() const override
  { return xs_.DiffusionInitialized(); }

  const bool ScatteringInitialized() const override
  { return xs_.ScatteringInitialized(); }

  const std::vector<double>& SigmaTotal() const override
  { return xs_.SigmaTotal(); }

  const std::vector<double>& SigmaAbsorption() const override
  { return xs_.SigmaAbsorption(); }

  const std::vector<double>& SigmaFission() const override
  { return xs_.SigmaFission(); }

  const std::vector<double>& NuSigmaF() const override
  { return xs_.NuSigmaF(); }

  const std::vector<double>& NuPromptSigmaF() const override
  { return xs_.NuPromptSigmaF(); }

  const std::vector<double>& NuDelayedSigmaF() const override
  { return xs_.NuDelayedSigmaF(); }

  const std::vector<double>& InverseVelocity() const override
  { return xs_.InverseVelocity(); }

  const std::vector<chi_math::SparseMatrix>& TransferMatrices() const override
  { return xs_.TransferMatrices(); }

  const chi_math::SparseMatrix& TransferMatrix(unsigned int ell) const override
  { return xs_.TransferMatrix(ell); }

  const std::vector<std::vector<double>> ProductionMatrix() const override
  { return xs_.ProductionMatrix(); }

  const std::vector<Precursor>& Precursors() const override
  { return xs_.Precursors(); }

  const std::vector<double>& DiffusionCoefficient() const override
  { return xs_.DiffusionCoefficient(); }

  const std::vector<double>& SigmaRemoval() const override
  { return xs_.SigmaRemoval(); }

  const std::vector<double>& SigmaSGtoG() const override
  { return xs_.SigmaSGtoG(); }

};


}//namespace chi_physics


#endif