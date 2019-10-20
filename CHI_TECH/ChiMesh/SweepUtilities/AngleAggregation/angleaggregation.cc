#include "angleaggregation.h"


//###################################################################
/** Gets the L^infinity norm of the relative change of the
 * delayed Psi values across either
 * intra-location or inter-location cyclic interfaces. */
double chi_mesh::sweep_management::AngleAggregation::GetDelayedPsiNorm()
{
  double loc_ret_val = 0.0;

  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      for (auto prelocI_norm : angset->delayed_prelocI_norm)
        loc_ret_val = std::max(prelocI_norm,loc_ret_val);

  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      loc_ret_val = std::max(angset->delayed_local_norm,loc_ret_val);

  double ret_val = 0.0;

  MPI_Allreduce(&loc_ret_val,&ret_val,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  return ret_val;
}

//###################################################################
/** Resets all the intra-location and inter-location cyclic interfaces.*/
void chi_mesh::sweep_management::AngleAggregation::ResetDelayedPsi()
{
  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      for (auto& delayed_data : angset->delayed_prelocI_outgoing_psi)
        delayed_data.assign(delayed_data.size(),0.0);

  for (auto angsetgrp : angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      angset->delayed_local_psi.assign(angset->delayed_local_psi.size(),0.0);
}