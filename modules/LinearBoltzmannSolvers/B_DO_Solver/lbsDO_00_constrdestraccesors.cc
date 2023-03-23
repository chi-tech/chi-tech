#include "lbs_discrete_ordinates_solver.h"

lbs::LBSDiscreteOrdinatesSolver::LBSDiscreteOrdinatesSolver(const std::string &text_name) :
  LBSSolver(text_name)
{
}

/**Destructor for LBS*/
lbs::LBSDiscreteOrdinatesSolver::~LBSDiscreteOrdinatesSolver()
{
  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }
}

/**Gets the local and global number of iterative unknowns. This normally is
 * only the flux moments, however, the sweep based solvers might include
 * delayed angular fluxes in this number.*/
std::pair<size_t, size_t> lbs::LBSDiscreteOrdinatesSolver::
GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_globl_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info =
      groupset.angle_agg_.GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_globl_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_globl_dofs = num_globl_phi_dofs + num_globl_psi_dofs;

  return {num_local_dofs, num_globl_dofs};
}