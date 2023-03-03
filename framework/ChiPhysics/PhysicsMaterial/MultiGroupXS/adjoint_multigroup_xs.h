#ifndef ADJOINT_MULTIGROUP_XS_H
#define ADJOINT_MULTIGROUP_XS_H

#include "multigroup_xs.h"


namespace chi_physics
{

/**
 * A class for accessing multi-group cross sections for adjoint problems.
 */
class AdjointMultiGroupXS : public MultiGroupXSBase
{
private:
  const MultiGroupXS& xs_;
  std::vector<chi_math::SparseMatrix> transposed_transfer_matrices_;
  std::vector<std::vector<double>> transposed_production_matrices_;

public:
  AdjointMultiGroupXS() = delete;
  AdjointMultiGroupXS(const AdjointMultiGroupXS&) = delete;
  AdjointMultiGroupXS(AdjointMultiGroupXS&&) = delete;

  explicit AdjointMultiGroupXS(const MultiGroupXS& xs);

  //Accessors
  const unsigned int NumGroups() const override { return xs_.NumGroups(); }

  const unsigned int ScatteringOrder() const override
  { return xs_.ScatteringOrder(); }

  const unsigned int NumPrecursors() const override
  { return xs_.NumPrecursors(); }

  const bool IsFissionable() const override { return xs_.IsFissionable(); }

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
  { return transposed_transfer_matrices_; }

  const chi_math::SparseMatrix& TransferMatrix(unsigned int ell) const override
  { return transposed_transfer_matrices_.at(ell); }

  const std::vector<std::vector<double>> ProductionMatrix() const override
  { return transposed_production_matrices_; }

  const std::vector<Precursor>& Precursors() const override
  { return xs_.Precursors(); }

  const std::vector<double>& DiffusionCoefficient() const override
  { return xs_.DiffusionCoefficient(); }

  const std::vector<double>& SigmaRemoval() const override
  { return xs_.SigmaRemoval(); }

  const std::vector<double>& SigmaSGtoG() const override
  { return xs_.SigmaSGtoG(); }

};

}

#endif //ADJOINT_MULTIGROUP_XS_H
