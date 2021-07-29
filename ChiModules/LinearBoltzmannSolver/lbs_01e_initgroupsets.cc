#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/**Initializes common groupset items.*/
void LinearBoltzmann::Solver::InitializeGroupsets()
{
  for (auto& groupset : groupsets)
  {
    //================================================== Build groupset angular
    //                                                   flux unknown manager
    groupset.psi_uk_man.unknowns.clear();
    size_t num_angles = groupset.quadrature->abscissae.size();
    size_t num_groups = groupset.groups.size();
    auto& grpset_psi_uk_man = groupset.psi_uk_man;

    for (unsigned int n=0; n<num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);

    //================================================== Setup group angular flux
    if (options.save_angular_flux)
    {
      size_t num_ang_unknowns = discretization->GetNumLocalDOFs(grpset_psi_uk_man);
      groupset.psi_to_be_saved = true;
      groupset.num_psi_unknowns_local = num_ang_unknowns;
      groupset.psi_new_local.assign(num_ang_unknowns,0.0);
    }

    groupset.BuildDiscMomOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildMomDiscOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildSubsets();
  }//for groupset
}