#include "angleaggregation.h"

#include "chi_log.h"
;

#include "chi_mpi.h"


//###################################################################
/** Sets up the angle-aggregation object. */
void chi_mesh::sweep_management::AngleAggregation::
  Setup(const std::vector<std::shared_ptr<SweepBndry>>& in_sim_boundaries,
        size_t in_number_of_groups,
        size_t in_number_of_group_subsets,
        std::shared_ptr<chi_math::AngularQuadrature> &in_quadrature,
        chi_mesh::MeshContinuumPtr& in_grid)
{
  sim_boundaries = in_sim_boundaries;
  number_of_groups = in_number_of_groups;
  number_of_group_subsets = in_number_of_group_subsets;
  quadrature = in_quadrature;
  grid = in_grid;

  is_setup = true;
}

//###################################################################
/** Resets all the outgoing intra-location and inter-location
 * cyclic interfaces.*/
void chi_mesh::sweep_management::AngleAggregation::ZeroOutgoingDelayedPsi()
{
  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.angle_sets)
      for (auto& delayed_data : angset->delayed_prelocI_outgoing_psi)
        delayed_data.assign(delayed_data.size(),0.0);

  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.angle_sets)
      angset->delayed_local_psi.assign(angset->delayed_local_psi.size(),0.0);
}

//###################################################################
/** Resets all the incoming intra-location and inter-location
 * cyclic interfaces.*/
void chi_mesh::sweep_management::AngleAggregation::ZeroIncomingDelayedPsi()
{
  //======================================== Opposing reflecting bndries
  for (const auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux_old)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = 0.0;

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& val : angle_set->delayed_local_psi_old)
        val = 0.0;

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi_old)
        for (auto& val : loc_vector)
          val = 0.0;
}

//###################################################################
/** Initializes reflecting boundary conditions. */
void chi_mesh::sweep_management::AngleAggregation::InitializeReflectingBCs()
{
  const double epsilon = 1.0e-8;

  int total_reflect_cells = 0;
  int total_reflect_faces = 0;
  int total_reflect_size = 0;

  bool reflecting_bcs_initialized=false;

  int bndry_id=0;
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      size_t tot_num_angles = quadrature->abscissae.size();
      size_t num_local_cells = grid->local_cell_glob_indices.size();
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      rbndry.reflected_anglenum.resize(tot_num_angles,-1);
      rbndry.angle_readyflags.resize(tot_num_angles,
                        std::vector<bool>(number_of_group_subsets,false));

      //========================================= Determine reflected angle
      for (int n=0; n<tot_num_angles; ++n)
      {
        auto omega_reflected = (quadrature->omegas[n]) -
                               rbndry.normal*
                               quadrature->omegas[n].Dot(rbndry.normal)*2.0;
        for (int nstar=0; nstar<tot_num_angles; ++nstar)
          if (omega_reflected.Dot(quadrature->omegas[nstar])> (1.0-epsilon))
            {rbndry.reflected_anglenum[n] = nstar;break;}

        if (rbndry.reflected_anglenum[n]<0)
        {
          chi::log.LogAllError()
            << "Reflected angle not found for angle " << n
            << " with direction " << quadrature->omegas[n].PrintS()
            << ". This can happen for two reasons: i) A quadrature is used"
               " that is not symmetric about the axis associated with the "
               "reflected boundary, or ii) the reflecting boundary is not "
               "aligned with any reflecting axis of the quadrature.";
          chi::Exit(EXIT_FAILURE);
        }
      }

      //========================================= For angles
      rbndry.hetero_boundary_flux.clear();
      rbndry.hetero_boundary_flux_old.clear();
      rbndry.hetero_boundary_flux.resize(tot_num_angles);
      for (int n=0; n<tot_num_angles; ++n)
      {
        //Only continue if omega is outgoing
        if ( quadrature->omegas[n].Dot(rbndry.normal)< 0.0 )
          continue;

        //================================== For cells
        auto& cell_vec = rbndry.hetero_boundary_flux[n];
        cell_vec.resize(num_local_cells);
        for (const auto& cell : grid->local_cells)
        {
          int c = cell.local_id;

          //=========================== Check cell on ref bndry
          bool on_ref_bndry = false;
          for (const auto& face : cell.faces){
            if ( (not face.has_neighbor) and
                 (face.normal.Dot(rbndry.normal) > 0.999999) )
            {
              on_ref_bndry = true;
              break;
            }
          }
          if (not on_ref_bndry) continue;
          total_reflect_cells += 1;

          //=========================== If cell on ref bndry
          cell_vec[c].resize(cell.faces.size());
          int f=0;
          for (const auto& face : cell.faces)
          {
            if ( (not face.has_neighbor) and
                 (face.normal.Dot(rbndry.normal) > 0.999999) )
            {
              cell_vec[c][f].clear();
              cell_vec[c][f].resize(face.vertex_ids.size(),
                                    std::vector<double>(number_of_groups,0.0));
              total_reflect_faces += 1;
              total_reflect_size  += face.vertex_ids.size()*number_of_groups;
            }
            ++f;
          }
        }//for cells
      }//for angles

      //========================================= Determine if boundary is
      //                                          opposing reflecting
      if ((bndry_id == 1) and (sim_boundaries[0]->IsReflecting()))
        rbndry.opposing_reflected = true;
      if ((bndry_id == 3) and (sim_boundaries[2]->IsReflecting()))
        rbndry.opposing_reflected = true;
      if ((bndry_id == 5) and (sim_boundaries[4]->IsReflecting()))
        rbndry.opposing_reflected = true;

      if (rbndry.opposing_reflected)
        rbndry.hetero_boundary_flux_old = rbndry.hetero_boundary_flux;

      reflecting_bcs_initialized = true;
    }//if reflecting

    ++bndry_id;
  }//for bndry

  if (reflecting_bcs_initialized)
    chi::log.Log0Verbose1() << "Reflecting boundary conditions initialized.";

}

//###################################################################
/** Returns a pair of numbers containing the number of
 * delayed angular unknowns both locally and globally, respectively. */
std::pair<size_t,size_t> chi_mesh::sweep_management::AngleAggregation::
  GetNumDelayedAngularDOFs()
{
  //======================================== Check if this is already developed
  if (num_ang_unknowns_avail)
    return number_angular_unknowns;

  //======================================== If not developed
  size_t local_ang_unknowns = 0;

  //======================================== Opposing reflecting bndries
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                local_ang_unknowns += dofvec.size();

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      local_ang_unknowns += angle_set->delayed_local_psi.size();

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi)
        local_ang_unknowns += loc_vector.size();



  size_t global_ang_unknowns = 0;
  MPI_Allreduce(&local_ang_unknowns,
                &global_ang_unknowns,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);

  number_angular_unknowns = {local_ang_unknowns,global_ang_unknowns};

  num_ang_unknowns_avail = true;
  return number_angular_unknowns;
}

//###################################################################
/** Assembles angular unknowns into the reference vector. */
void chi_mesh::sweep_management::AngleAggregation::
  AppendDelayedAngularDOFsToArray(int &index, double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                {index++; x_ref[index] = val;}

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto val : angle_set->delayed_local_psi)
      {index++; x_ref[index] = val;}

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi)
        for (auto val : loc_vector)
        {index++; x_ref[index] = val;}
}

//###################################################################
/** Assembles angular unknowns into the reference vector. */
void chi_mesh::sweep_management::AngleAggregation::
  SetDelayedAngularDOFsFromArray(int &index, const double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux_old)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                {index++; val = x_ref[index];}

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& val : angle_set->delayed_local_psi_old)
      {index++; val = x_ref[index];}

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi_old)
        for (auto& val : loc_vector)
        {index++; val = x_ref[index];}
}

//###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
std::vector<double> chi_mesh::sweep_management::AngleAggregation::
  GetDelayedAngularDOFsAsSTLVector()
{
  std::vector<double> psi_vector;

  auto psi_size = GetNumDelayedAngularDOFs();
  psi_vector.reserve(psi_size.first);

  //======================================== Opposing reflecting bndries
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                  psi_vector.push_back(val);

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto val : angle_set->delayed_local_psi)
        psi_vector.push_back(val);

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi)
        for (auto val : loc_vector)
          psi_vector.push_back(val);

  return psi_vector;
}

//###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
void chi_mesh::sweep_management::AngleAggregation::
  SetDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector)
{
  auto psi_size = GetNumDelayedAngularDOFs();
  size_t stl_size = stl_vector.size();
  if (stl_size != psi_size.first)
    throw std::logic_error(std::string(__FUNCTION__) + ": STL-vector size "
                           "is incompatible with number angular unknowns stored "
                           "in the angle-aggregation object.");

  size_t index = 0;
  //======================================== Opposing reflecting bndries
  for (auto& bndry : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.opposing_reflected)
        for (auto& angle : rbndry.hetero_boundary_flux_old)
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = stl_vector[index++];

    }//if reflecting
  }//for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& val : angle_set->delayed_local_psi_old)
        val = stl_vector[index++];

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.angle_sets)
      for (auto& loc_vector : angle_set->delayed_prelocI_outgoing_psi_old)
        for (auto& val : loc_vector)
          val = stl_vector[index++];
}