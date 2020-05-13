#include "lbs_linear_boltzman_solver.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Performs general input checks before initialization continues.*/
void LinearBoltzman::Solver::PerformInputChecks()
{
  if (groups.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzman::Solver: No groups added to solver.";
    exit(EXIT_FAILURE);
  }
  if (group_sets.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzman::Solver: No group-sets added to solver.";
    exit(EXIT_FAILURE);
  }
  int grpset_counter=0;
  for (auto group_set : group_sets)
  {
    if (group_set->groups.empty())
    {
      chi_log.Log(LOG_ALLERROR)
        << "LinearBoltzman::Solver: No groups added to groupset "
        << grpset_counter << ".";
      exit(EXIT_FAILURE);
    }
    ++grpset_counter;
  }
  if (discretization == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzman::Solver: No discretization method set.";
    exit(EXIT_FAILURE);
  }
  if (regions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzman::Solver: No regions added to solver.";
    exit(EXIT_FAILURE);
  }
  chi_mesh::Region*  aregion = regions.back();
  grid                       = aregion->GetGrid();


}