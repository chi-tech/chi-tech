#include "lbs_sweepchunk.h"

namespace lbs
{
const double* SweepSurfaceStatusInfo::GetUpwindPsi(int fj) const
{
  const double* psi;
  if (on_local_face)
    psi = fluds->UpwindPsi(spls_index, in_face_counter, fj, 0, angle_set_index);
  else if (not on_boundary)
    psi = fluds->NLUpwindPsi(preloc_face_counter, fj, 0, angle_set_index);
  else
    psi = angle_set->PsiBndry(bndry_id,
                              angle_num,
                              cell_local_id,
                              f,
                              fj,
                              gs_gi_,
                              gs_ss_begin_,
                              surface_source_active);
  return psi;
}

double* SweepSurfaceStatusInfo::GetDownwindPsi(int fi) const
{
  double* psi;
  if (on_local_face)
    psi = fluds->OutgoingPsi(spls_index, out_face_counter, fi, angle_set_index);
  else if (not on_boundary)
    psi = fluds->NLOutgoingPsi(deploc_face_counter, fi, angle_set_index);
  else if (is_reflecting_bndry_)
    psi = angle_set->ReflectingPsiOutBoundBndry(
      bndry_id, angle_num, cell_local_id, f, fi, gs_ss_begin_);
  else
    psi = nullptr;

  return psi;
}
} // namespace lbs