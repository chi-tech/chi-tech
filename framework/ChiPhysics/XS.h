#ifndef CHITECH_XS_H
#define CHITECH_XS_H

#include <vector>
#include <cstddef>
#include <memory>

namespace chi_math
{
  class SparseMatrix;
}

namespace chi_physics
{

struct XSData
{
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
};

class XSBase
{
private:
  std::unique_ptr<XSData> data = nullptr;
  bool initialized = false;

public:
  virtual const std::vector<double>& Sigma_t() const = 0;
  virtual const std::vector<double>& Sigma_capture() const = 0;
  virtual const std::vector<double>& Sigma_f() const = 0;
  virtual const std::vector<double>& Nu() const = 0;
  //etc.

  virtual const std::vector<chi_math::SparseMatrix>& Sigma_s_Matrices() const = 0;
  virtual const chi_math::SparseMatrix& Sigma_f_Matrix() const = 0;
  //etc.
};

/**
 * Other stuff...
 *
 * Acronyms and definitions:
 * dofs = Degress of freedom
 * exclusive = no ghosts
 * inclusive = includes ghosts
 * */
struct GenericDoFsPerObject
{
  std::vector<std::vector<unsigned int>> dofs_per_object_exclusive;
  std::vector<std::vector<unsigned int>> dofs_per_object_inclusive;
  std::vector<std::vector<unsigned int>> object_index;
  std::vector<std::vector<unsigned int>> first_object_index_on_face;
};

}//namespace chi_physics

#endif //CHITECH_XS_H
