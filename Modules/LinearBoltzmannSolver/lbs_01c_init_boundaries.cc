#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/**Initializes transport related boundaries. */
void LinearBoltzmann::Solver::InitializeBoundaries()
{
  //================================================== Initialize default
  //                                                   incident boundary
  typedef chi_mesh::sweep_management::BoundaryVacuum SweepVacuumBndry;
  typedef chi_mesh::sweep_management::BoundaryIncidentHomogenous SweepIncHomoBndry;
  typedef chi_mesh::sweep_management::BoundaryReflecting SweepReflectingBndry;
  std::vector<std::vector<double>>& flux_vec = incident_P0_mg_boundaries;

  // Defining default Vacuum boundary
  zero_boundary.resize(groups.size(),0.0);

  // ================================================= Populate boundaries
  if (sweep_boundaries.empty())
  {
    chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
    chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
    chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
    int bndry_id=0;
    for (auto bndry_type : boundary_types)
    {
      int vec_index = bndry_type.second;

      if (bndry_type.first == LinearBoltzmann::BoundaryType::VACUUM)
        sweep_boundaries.push_back(new SweepVacuumBndry(zero_boundary));
      else if (bndry_type.first == LinearBoltzmann::BoundaryType::INCIDENT_ISOTROPIC)
        sweep_boundaries.push_back(new SweepIncHomoBndry(flux_vec[vec_index]));
      else if (bndry_type.first == LinearBoltzmann::BoundaryType::REFLECTING)
      {
        chi_mesh::Normal normal;
        if (bndry_id == 0) normal = ihat;
        if (bndry_id == 1) normal = ihat*-1.0;
        if (bndry_id == 2) normal = jhat;
        if (bndry_id == 3) normal = jhat*-1.0;
        if (bndry_id == 4) normal = khat;
        if (bndry_id == 5) normal = khat*-1.0;

        sweep_boundaries.push_back(
          new SweepReflectingBndry(zero_boundary, normal));
      }

      ++bndry_id;
    }
  }//if empty

}