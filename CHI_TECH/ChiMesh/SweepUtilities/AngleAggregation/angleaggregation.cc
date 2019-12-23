#include "angleaggregation.h"

#include "chi_log.h"
extern ChiLog chi_log;

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

//###################################################################
/** Initializes reflecting boundary conditions. */
void chi_mesh::sweep_management::AngleAggregation::InitializeReflectingBCs()
{
  chi_log.Log(LOG_0) << "Initializing Reflecting boundary conditions.";
  int total_reflect_cells = 0;
  int total_reflect_faces = 0;
  int total_reflect_size = 0;
  for (auto bndry : sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      size_t tot_num_angles = quadrature->abscissae.size();
      size_t num_local_cells = grid->local_cell_glob_indices.size();
      auto rbndry = (chi_mesh::sweep_management::BoundaryReflecting*)bndry;

      rbndry->reflected_anglenum.resize(tot_num_angles,-1);
      rbndry->angle_readyflags.resize(tot_num_angles,
                        std::vector<bool>(number_of_group_subsets,false));

      //========================================= For angles
      rbndry->hetero_boundary_flux.resize(tot_num_angles);
      for (int n=0; n<tot_num_angles; ++n)
      {
        //================================== Determine reflected angles
        //Only continue if omega is outgoing
        if ( quadrature->omegas[n]->Dot(rbndry->normal)> 0.0 )
        {
          auto omega_reflected = (*quadrature->omegas[n]) -
                                 rbndry->normal*
                                 quadrature->omegas[n]->Dot(rbndry->normal)*2.0;
          for (int nstar=0; nstar<tot_num_angles; ++nstar)
            if (omega_reflected.Dot(*quadrature->omegas[nstar])> 0.999999){
              rbndry->reflected_anglenum[n] = nstar;break;}
        } else
        {
          auto omega_reflected = (*quadrature->omegas[n]) -
                                 rbndry->normal*
                                 quadrature->omegas[n]->Dot(rbndry->normal)*2.0;
          for (int nstar=0; nstar<tot_num_angles; ++nstar)
            if (omega_reflected.Dot(*quadrature->omegas[nstar])> 0.999999){
              rbndry->reflected_anglenum[n] = nstar;break;}
          continue;
        }

        //================================== For cells
        auto& cell_vec = rbndry->hetero_boundary_flux[n];
        cell_vec.resize(num_local_cells);
        for (size_t c=0; c<num_local_cells; ++c)
        {
          auto cell = grid->cells[grid->local_cell_glob_indices[c]];

          //=========================== Check cell on ref bndry
          bool on_ref_bndry = false;
          for (auto& face : cell->faces){
            if ( (face.neighbor < 0) and
                 (face.normal.Dot(rbndry->normal) > 0.999999) )
            {
              on_ref_bndry = true;
              break;
            }
          }
          if (not on_ref_bndry) continue;
          total_reflect_cells += 1;

          //=========================== If cell on ref bndry
          cell_vec[c].resize(cell->faces.size());
          int f=0;
          for (auto& face : cell->faces)
          {
            if ( (face.neighbor < 0) and
                 (face.normal.Dot(rbndry->normal) > 0.999999) )
            {
              cell_vec[c][f].resize(face.vertex_ids.size(),
                                    std::vector<double>(number_of_groups,0.0));
              total_reflect_faces += 1;
              total_reflect_size  += face.vertex_ids.size()*number_of_groups;
            }
            ++f;
          }
        }//for cells
      }//for angles
    }//if reflecting
  }//for bndry

  chi_log.Log(LOG_0) << "Total reflecting cells = " << total_reflect_cells;
  chi_log.Log(LOG_0) << "Total reflecting faces = " << total_reflect_faces;
  chi_log.Log(LOG_0) << "Total reflecting variables = " << total_reflect_size;
}

//###################################################################
/** Initializes reflecting boundary conditions. */
void chi_mesh::sweep_management::AngleAggregation::ResetReflectingBCs()
{
  for (auto bndry : sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      auto rbndry = (chi_mesh::sweep_management::BoundaryReflecting*)bndry;

      for (auto& angle : rbndry->hetero_boundary_flux)
        for (auto& cellvec : angle)
          for (auto& facevec : cellvec)
            for (auto& dofvec : facevec)
              for (auto& val : dofvec)
                val = 0.0;
    }//if reflecting
  }//for bndry
}