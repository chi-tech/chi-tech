#ifndef CHITECH_LBS_SOURCE_CONTEXT_H
#define CHITECH_LBS_SOURCE_CONTEXT_H

#include "ChiPhysics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"

#include <cstddef>
#include <vector>


namespace chi_math
{
  class SparseMatrix;
}

namespace lbs
{
class LBSSolver;
class LBSGroupset;

class SourceContext
{
protected:
  bool apply_wgs_scatter_src_ = false;
  bool apply_ags_scatter_src_ = false;
  bool apply_wgs_fission_src_ = false;
  bool apply_ags_fission_src_ = false;
  bool suppress_wg_scatter_src_ = false;

  size_t gs_i_      = 0;
  size_t gs_f_      = 0;
  size_t first_grp_ = 0;
  size_t last_grp_  = 0;

  double cell_volume_ = 0.0;
  size_t g_         = 0;

  const double* fixed_src_moments_ = nullptr;

public:
  SourceContext() = default;

  //=================================== Setters
  void SetFlags(bool apply_wgs_scatter_src,
                bool apply_ags_scatter_src,
                bool apply_wgs_fission_src,
                bool apply_ags_fission_src,
                bool suppress_wg_scatter_src);

  void SetCellVolume(double cell_volume);

  void SetGroupsetBounds(size_t gs_i, size_t gs_f);

  void SetGroupBounds(size_t first_grp, size_t last_grp);

  void SetGroupIndex(size_t g);

  void SetFixedSrcMomentsData(const double* fixed_src_moments);

  //=================================== Getters
  std::pair<size_t, size_t> GetGroupsetBounds() const;

  std::pair<size_t, size_t> GetGroupBounds() const;

  size_t GetActiveGroupIndex() const;

  const double* GetFixedSrcMomentsData() const;

  //=================================== Virtual Methods
  virtual double AddSourceMoments() const = 0;

  virtual
  double AddScattering(const chi_math::SparseMatrix& S,
                       const double* phi) const = 0;

  virtual
  double AddPromptFission(const std::vector<double>& F_g,
                          const double* phi) const = 0;

  typedef std::vector<chi_physics::MultiGroupXS::Precursor> PrecursorList;
  virtual
  double AddDelayedFission(const PrecursorList& precursors,
                           const std::vector<double>& nu_delayed_sigma_f,
                           const double* phi) const = 0;
};

}//namespace lbs

#endif //CHITECH_LBS_SOURCE_CONTEXT_H
