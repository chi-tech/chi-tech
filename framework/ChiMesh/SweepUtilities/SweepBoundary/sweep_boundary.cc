#include "sweep_boundaries.h"

#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double* chi_mesh::sweep_management::SweepBoundary::
HeterogeneousPsiIncoming(uint64_t cell_local_id,
                         int face_num,
                         int fi,
                         int angle_num,
                         int group_num,
                         int gs_ss_begin)
{
  chi::log.LogAllError()
    << "HeterogeneousPsiIncoming call made to boundary "
       "that has no such information.";
  chi::Exit(EXIT_FAILURE);
  return nullptr;
}

//###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double* chi_mesh::sweep_management::SweepBoundary::
HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                         int face_num,
                         int fi,
                         int angle_num,
                         int gs_ss_begin)
{
  chi::log.LogAllError()
    << "HeterogeneousPsiOutgoing call made to boundary "
       "that has no such information.";
  chi::Exit(EXIT_FAILURE);
  return nullptr;
}