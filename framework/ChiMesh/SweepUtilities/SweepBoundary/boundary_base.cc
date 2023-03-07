#include "sweep_boundaries.h"

#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryBase::
HeterogenousPsiIncoming(uint64_t cell_local_id,
                        int face_num,
                        int fi,
                        int angle_num,
                        int group_num,
                        int gs_ss_begin)
{
  chi::log.LogAllError()
    << "HeterogenousPsiIncoming call made to boundary "
       "that has no such information.";
  chi::Exit(EXIT_FAILURE);
  return nullptr;
}

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryBase::
HeterogenousPsiOutgoing(uint64_t cell_local_id,
                        int face_num,
                        int fi,
                        int angle_num,
                        int gs_ss_begin)
{
  chi::log.LogAllError()
    << "HeterogenousPsiOutgoing call made to boundary "
       "that has no such information.";
  chi::Exit(EXIT_FAILURE);
  return nullptr;
}