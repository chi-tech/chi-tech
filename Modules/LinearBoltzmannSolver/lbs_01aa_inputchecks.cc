#include "lbs_linear_boltzmann_solver.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Performs general input checks before initialization continues.*/
void LinearBoltzmann::Solver::PerformInputChecks()
{
  if (groups.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzmann::Solver: No groups added to solver.";
    exit(EXIT_FAILURE);
  }
  if (group_sets.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzmann::Solver: No group-sets added to solver.";
    exit(EXIT_FAILURE);
  }
  int grpset_counter=0;
  for (auto& group_set : group_sets)
  {
    if (group_set.groups.empty())
    {
      chi_log.Log(LOG_ALLERROR)
        << "LinearBoltzmann::Solver: No groups added to groupset "
        << grpset_counter << ".";
      exit(EXIT_FAILURE);
    }
    ++grpset_counter;
  }
  if (discretization == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzmann::Solver: No discretization method set.";
    exit(EXIT_FAILURE);
  }
  if (regions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzmann::Solver: No regions added to solver.";
    exit(EXIT_FAILURE);
  }
  chi_mesh::Region*  aregion = regions.back();
  grid                       = aregion->GetGrid();

  if (grid == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "LinearBoltzmann::Solver: No grid available from region.";
    exit(EXIT_FAILURE);
  }

  //======================================== Determine geometry type
  if (grid->local_cells[0].Type() == chi_mesh::CellType::SLAB)
  {
    options.geometry_type = GeometryType::ONED_SLAB;
  }
  else if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYGON)
  {
    options.geometry_type = GeometryType::TWOD_CARTESIAN;
  }
  else if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYHEDRON)
  {
    options.geometry_type = GeometryType::THREED_CARTESIAN;
  }
}