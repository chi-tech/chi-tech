#include "lbs_linear_boltzmann_solver.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_mpi.h>
#include <chi_log.h>


extern ChiLog& chi_log;

#include <iomanip>
#include "ChiConsole/chi_console.h"

extern ChiConsole&  chi_console;



//###################################################################
/** Initialize the solver.*/
void lbs::SteadySolver::Initialize()
{
  PerformInputChecks(); //a
  PrintSimHeader(); //b
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Add unique material ids
  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (auto& cell : grid->local_cells)
  {
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id<0)
      ++invalid_mat_cell_count;
  }

  if (invalid_mat_cell_count>0)
  {
    chi_log.Log(LOG_ALLWARNING)
      << "Number of invalid material cells: " << invalid_mat_cell_count;
  }

  //================================================== Initialize materials
  InitMaterials(unique_material_ids); //c

  //================================================== Init spatial discretization
  InitializeSpatialDiscretization(); //d

  //================================================== Initialize groupsets
  InitializeGroupsets(); //e

  //================================================== Compute n. moments
  ComputeNumberOfMoments(); //f

  //================================================== Initialize parrays
  chi_log.Log(LOG_0)
    << "Initializing parallel arrays. " << std::endl;

  InitializeParrays();//g

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Done with parallel arrays.                Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB" << std::endl;

  //================================================== Initialize boundaries
  InitializeBoundaries();//h

  //================================================== Initialize sources
  InitializePointSources();

}
