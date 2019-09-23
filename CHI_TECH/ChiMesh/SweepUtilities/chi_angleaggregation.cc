#include "chi_angleaggregation.h"



double chi_mesh::SweepManagement::AngleAggregation::GetDelayedPsiNorm()
{
  double loc_ret_val = 0.0;

  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      for (auto prelocI_norm : angset->delayed_prelocI_norm)
        loc_ret_val += prelocI_norm;

  double ret_val = 0.0;

  MPI_Allreduce(&loc_ret_val,&ret_val,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return ret_val;
}

void chi_mesh::SweepManagement::AngleAggregation::ResetDelayedPsi()
{
  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      for (auto delayed_data : angset->delayed_prelocI_outgoing_psi)
        delayed_data.assign(delayed_data.size(),0.0);
}