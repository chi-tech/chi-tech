#include "sweep_boundaries.h"

//###################################################################
/**Returns a pointer to a homogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryIsotropicHomogenous::
HeterogeneousPsiIncoming(uint64_t cell_local_id,
                         int face_num,
                         int fi,
                         int angle_num,
                         int group_num,
                         int gs_ss_begin)
{
  return &boundary_flux[group_num];
}