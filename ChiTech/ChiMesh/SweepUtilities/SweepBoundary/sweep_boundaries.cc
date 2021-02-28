#include "sweep_boundaries.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;

////###################################################################
///**Returns a flag indicating whether this bndry is reflecting or not.*/
//bool chi_mesh::sweep_management::BoundaryBase::IsReflecting()
//{
//  if (this->Type() == BoundaryType::REFLECTING)
//    return true;
//  else
//    return false;
//}

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryBase::
HeterogenousPsiIncoming(
                int angle_num,
                int cell_local_id,
                int face_num,
                int fi,
                int gs_ss_begin)
{
  chi_log.Log(LOG_ALLERROR)
    << "HeterogenousPsiIncoming call made to boundary "
       "that has no such information.";
  exit(EXIT_FAILURE);
}

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryBase::
HeterogenousPsiOutgoing(
                int angle_num,
                int cell_local_id,
                int face_num,
                int fi,
                int gs_ss_begin)
{
  chi_log.Log(LOG_ALLERROR)
    << "HeterogenousPsiOutgoing call made to boundary "
       "that has no such information.";
  exit(EXIT_FAILURE);
}

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryReflecting::
HeterogenousPsiIncoming(
                int angle_num,
                int cell_local_id,
                int face_num,
                int fi,
                int gs_ss_begin)
{
  double* Psi = zero_boundary_flux.data();

  int reflected_angle_num = reflected_anglenum[angle_num];

  if (opposing_reflected)
  {
    Psi = &hetero_boundary_flux_old[reflected_angle_num]
                                   [cell_local_id]
                                   [face_num]
                                   [fi][gs_ss_begin];
  }
  else
  {
    Psi = &hetero_boundary_flux[reflected_angle_num]
                               [cell_local_id]
                               [face_num]
                               [fi][gs_ss_begin];
  }

  return Psi;
}

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryReflecting::
HeterogenousPsiOutgoing(
  int angle_num,
  int cell_local_id,
  int face_num,
  int fi,
  int gs_ss_begin)
{
  return &hetero_boundary_flux[angle_num]
                              [cell_local_id]
                              [face_num]
                              [fi][gs_ss_begin];;
}


//###################################################################
/**Sets flags indicating reflected angles are ready to execute.*/
void chi_mesh::sweep_management::BoundaryReflecting::
  UpdateAnglesReadyStatus(std::vector<int> angles, int gs_ss)
{
  for (auto& n : angles)
    angle_readyflags[reflected_anglenum[n]][gs_ss] = true;
}

//###################################################################
/**Checks to see if angles are ready to execute.*/
bool chi_mesh::sweep_management::BoundaryReflecting::
  CheckAnglesReadyStatus(std::vector<int> angles, int gs_ss)
{
  if (opposing_reflected) return true;
  bool ready_flag = true;
  for (auto& n : angles)
    if (hetero_boundary_flux[reflected_anglenum[n]].size()>0)
      if (not angle_readyflags[n][gs_ss]) return false;

  return ready_flag;
}

//###################################################################
/**Resets angle ready flags to false.*/
void chi_mesh::sweep_management::BoundaryReflecting::
  ResetAnglesReadyStatus()
{
  double local_pw_change = 0.0;
  if (opposing_reflected)
  {
    int n=0;
    for (auto& angle : hetero_boundary_flux){
      int c=0;
      for (auto& cellvec : angle){
        int f=0;
        for (auto& facevec : cellvec){
          int i=0;
          for (auto& dofvec : facevec){
            int g=0;
            for (auto& new_val : dofvec)
            {
              double old_val = hetero_boundary_flux_old[n][c][f][i][g];
              double delta_val = std::fabs(new_val-old_val);
              double max_val = std::max(new_val,old_val);

              if (max_val >= std::numeric_limits<double>::min())
                local_pw_change = std::max(delta_val/max_val,local_pw_change);
              else
                local_pw_change = std::max(delta_val,local_pw_change);

              hetero_boundary_flux_old[n][c][f][i][g] = new_val;
              ++g;
            }
            ++i;
          }
          ++f;
        }
        ++c;
      }
      ++n;
    }
    MPI_Allreduce(&local_pw_change,&pw_change,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }

  for (auto& flags : angle_readyflags)
    for (int gs_ss=0; gs_ss<flags.size(); ++gs_ss)
      flags[gs_ss] = false;
}