#include "lbs_linear_boltzmann_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

;

//###################################################################
/**Performs general input checks before initialization continues.*/
void lbs::SteadySolver::PerformInputChecks()
{
  if (groups.empty())
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No groups added to solver.";
    exit(EXIT_FAILURE);
  }
  if (groupsets.empty())
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No group-sets added to solver.";
    exit(EXIT_FAILURE);
  }
  int grpset_counter=0;
  for (auto& group_set : groupsets)
  {
    if (group_set.groups.empty())
    {
      chi::log.LogAllError()
        << "LinearBoltzmann::Solver: No groups added to groupset "
        << grpset_counter << ".";
      exit(EXIT_FAILURE);
    }
    ++grpset_counter;
  }
  if (options.sd_type == chi_math::SpatialDiscretizationType::UNDEFINED)
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No discretization method set.";
    exit(EXIT_FAILURE);
  }
//  if (regions.empty())
//  {
//    chi::log.LogAllError()
//      << "LinearBoltzmann::Solver: No regions added to solver.";
//    exit(EXIT_FAILURE);
//  }
//  chi_mesh::Region*  aregion = regions.back();
//  grid                       = aregion->GetGrid();

  grid = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid == nullptr)
  {
    chi::log.LogAllError()
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