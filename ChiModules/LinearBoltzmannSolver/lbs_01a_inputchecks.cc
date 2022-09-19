#include "lbs_linear_boltzmann_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

//###################################################################
/**Performs general input checks before initialization continues.*/
void lbs::SteadySolver::PerformInputChecks()
{
  if (groups.empty())
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No groups added to solver.";
    chi::Exit(EXIT_FAILURE);
  }
  if (groupsets.empty())
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No group-sets added to solver.";
    chi::Exit(EXIT_FAILURE);
  }
  int grpset_counter=0;
  for (auto& group_set : groupsets)
  {
    if (group_set.groups.empty())
    {
      chi::log.LogAllError()
        << "LinearBoltzmann::Solver: No groups added to groupset "
        << grpset_counter << ".";
      chi::Exit(EXIT_FAILURE);
    }
    ++grpset_counter;
  }
  if (options.sd_type == chi_math::SpatialDiscretizationType::UNDEFINED)
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No discretization method set.";
    chi::Exit(EXIT_FAILURE);
  }

  grid = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid == nullptr)
  {
    chi::log.LogAllError()
      << "LinearBoltzmann::Solver: No grid available from region.";
    chi::Exit(EXIT_FAILURE);
  }

  //======================================== Determine geometry type
  using namespace chi_mesh;
  const auto grid_attribs = grid->Attributes();
  if (grid_attribs & DIMENSION_1)
    options.geometry_type = GeometryType::ONED_SLAB;
  if (grid_attribs & DIMENSION_2)
    options.geometry_type = GeometryType::TWOD_CARTESIAN;
  if (grid_attribs & DIMENSION_3)
    options.geometry_type = GeometryType::THREED_CARTESIAN;

}