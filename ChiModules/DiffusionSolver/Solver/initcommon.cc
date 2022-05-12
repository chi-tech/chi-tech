#include "diffusion_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"


//###################################################################
/**Initialization of common to all solver types.*/
void chi_diffusion::Solver::InitializeCommonItems()
{
//  if (regions.empty())
//  {
//    chi_log.Log(LOG_ALLERROR)
//      << "chi_diffusion::Solver::InitializeCommonItems: No region added to solver.";
//    exit(EXIT_FAILURE);
//  }
//
//  auto& region = regions.back();
//  grid = region->GetGrid();

  grid = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " No grid defined.");

  auto globl_unique_bndry_ids = grid->GetDomainUniqueBoundaryIDs();

  uint64_t max_boundary_id = 0;
  for (const auto& id : globl_unique_bndry_ids)
    max_boundary_id = std::max(id,max_boundary_id);

  chi_log.Log() << "Max boundary id identified: " << max_boundary_id;

  for (int bndry=0; bndry<(max_boundary_id+1); bndry++)
  {
    if (boundary_preferences.find(bndry) != boundary_preferences.end())
    {
      BoundaryInfo bndry_info = boundary_preferences.at(bndry);
      auto& bndry_vals = bndry_info.second;
      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting:
        {
          boundaries.push_back(
            new chi_diffusion::BoundaryReflecting);
          chi_log.Log() << "Boundary " << bndry << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1,0.0);
          boundaries.push_back(
            new chi_diffusion::BoundaryDirichlet(bndry_vals[0]));
          chi_log.Log() << "Boundary " << bndry << " set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size()<3) bndry_vals.resize(3,0.0);
          boundaries.push_back(
            new chi_diffusion::BoundaryRobin(bndry_vals[0],
                                             bndry_vals[1],
                                             bndry_vals[2]));
          chi_log.Log() << "Boundary " << bndry << " set to robin.";
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries.push_back(new chi_diffusion::BoundaryRobin(0.25,0.5,0.0));
          chi_log.Log() << "Boundary " << bndry << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size()<3) bndry_vals.resize(3,0.0);
          boundaries.push_back(
            new chi_diffusion::BoundaryRobin(bndry_vals[0],
                                             bndry_vals[1],
                                             bndry_vals[2]));
          chi_log.Log() << "Boundary " << bndry << " set to neumann.";
          break;
        }
      }//switch boundary type
    }
    else
    {
      boundaries.push_back(new chi_diffusion::BoundaryDirichlet);
      chi_log.Log(LOG_0VERBOSE_1)
        << "No boundary preference found for boundary index " << bndry
        << "Dirichlet boundary added with zero boundary value.";
    }
  }//for bndry

  common_items_initialized = true;
}


