#include "../chi_mesh.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <algorithm>

//###################################################################
/** This method searches ref_glob_dep[ref_locI] for a dependency on
 * master_locI and if found will set found=true.
 * Its search is registered in search_history to prevent
 * cyclic search patterns. Optionally the routine can also be allowed to
 * search recursively.*/
void
chi_mesh::sweep_management::RecursivelyFindLocIDependency(
  std::vector<std::vector<int>>& ref_glob_dep,
  std::vector<int>& search_history,
  int master_locI,
  int ref_locI,
  bool& found,
  bool allow_recursive_search)
{
  if (found) return;
  if (ref_locI < 0) return;

  auto already_searched = (
    std::find(search_history.begin(),search_history.end(),ref_locI) !=
    search_history.end() );

  if (already_searched) return;
  search_history.push_back(ref_locI);


  for (auto dependency : ref_glob_dep[ref_locI])
  {
    if (dependency == master_locI)
    {
      found = true;
      break;
    }
    if (allow_recursive_search)
      RecursivelyFindLocIDependency(ref_glob_dep,
                                    search_history,
                                    master_locI,
                                    dependency,
                                    found,
                                    allow_recursive_search);
  }
}

//###################################################################
/** This method firstly loops over each location and then each location's
 * dependencies. It can optionally also go deeper and recursively loop.
 *
 * The main purpose is to find whether dependencies reference back to this
 * location, thereby identifying a cyclic dependency. When a cyclic
 * dependency is found it will be removed and the appropriate registers
 * will be set.*/
void
chi_mesh::sweep_management::RemoveGlobalCyclicDependencies(
  chi_mesh::sweep_management::SPDS *sweep_order,
  std::vector<std::vector<int>> &global_dependencies,
  bool allow_recursive_search,
  bool allow_cycles)
{
  std::vector<int> search_history;
  search_history.reserve(chi_mpi.location_id);
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {

    int c=-1;
    for (auto rlocI : global_dependencies[locI])
    {
      ++c;
      if (rlocI<0) continue;

      bool depends_on_locI = false;
      search_history.clear();
      search_history.reserve(chi_mpi.location_id);
      RecursivelyFindLocIDependency(global_dependencies,
                                    search_history,
                                    locI,
                                    rlocI,
                                    depends_on_locI,
                                    allow_recursive_search);

      if (depends_on_locI and allow_cycles)
      {
        global_dependencies[locI][c]  = -1;
        if (locI == chi_mpi.location_id)
        {
          auto dependent_location =
            std::find(sweep_order->location_dependencies.begin(),
                      sweep_order->location_dependencies.end(),
                      rlocI);
          sweep_order->location_dependencies.erase(dependent_location);
          sweep_order->delayed_location_dependencies.push_back(rlocI);
        }

        if (rlocI == chi_mpi.location_id)
        {
          sweep_order->delayed_location_successors.push_back(locI);
        }
      }
      else if (depends_on_locI and (not allow_cycles))
      {
        chi_log.Log(LOG_ALLERROR)
          << "Global cyclic dependency detected. This must be allowed"
             " by client applications.";
        exit(EXIT_FAILURE);
      }



    }//for locI dependency c
  }//for locI
}