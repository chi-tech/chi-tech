#include "source_context.h"

/**Sets current source scope flags.*/
void lbs::SourceContext::SetFlags(bool apply_wgs_scatter_src,
                                  bool apply_ags_scatter_src,
                                  bool apply_wgs_fission_src,
                                  bool apply_ags_fission_src,
                                  bool suppress_wg_scatter_src)
{
  apply_wgs_scatter_src_ = apply_wgs_scatter_src;
  apply_ags_scatter_src_ = apply_ags_scatter_src;
  apply_wgs_fission_src_ = apply_wgs_fission_src;
  apply_ags_fission_src_ = apply_ags_fission_src;
  suppress_wg_scatter_src_ = suppress_wg_scatter_src;
}

/**Sets a pair of numbers. The first number is the group index of where
 * the group index starts for the groupset, and the second number is the
 * group index of the last group in the groupset.*/
void lbs::SourceContext::SetGroupsetBounds(size_t gs_i, size_t gs_f)
{
  gs_i_ = gs_i;
  gs_f_ = gs_f;
}

/**Sets a pair of numbers. The first number is the group index of where
 * the group index starts for the problem, and the second number is the group
 * index of the last group in the problem.*/
void lbs::SourceContext::SetGroupBounds(size_t first_grp, size_t last_grp)
{
  first_grp_ = first_grp;
  last_grp_ = last_grp;
}

/**Set active cell's volume.*/
void lbs::SourceContext::SetCellVolume(double cell_volume)
{
  cell_volume_ = cell_volume;
}


/**Sets the currently active group.*/
void lbs::SourceContext::SetGroupIndex(size_t g)
{
  g_ = g;
}

/**Sets a pointer to the currently active read-only data.*/
void lbs::SourceContext::SetFixedSrcMomentsData(const double *fixed_src_moments)
{
  fixed_src_moments_ = fixed_src_moments;
}





/**Returns a pair of numbers. The first number is the group index of where
 * the group index starts for the groupset, and the second number is the
 * group index of the last group in the groupset.*/
std::pair<size_t, size_t> lbs::SourceContext::GetGroupsetBounds() const
{
  return {gs_i_, gs_f_};
}

/**Returns a pair of numbers. The first number is the group index of where
 * the group index starts for the problem, and the second number is the group
 * index of the last group in the problem.*/
std::pair<size_t, size_t> lbs::SourceContext::GetGroupBounds() const
{
  return {first_grp_, last_grp_};
}

/**Returns currently active group.*/
size_t lbs::SourceContext::GetActiveGroupIndex() const
{
  return g_;
}

/**Returns a pointer to the currently active read-only data.*/
const double* lbs::SourceContext::GetFixedSrcMomentsData() const
{
  return fixed_src_moments_;
}