#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_runtime.h"
#include "chi_log.h"
;

#include "ChiConsole/chi_console.h"



#include <iomanip>

void lbs::SteadySolver::InitializeSpatialDiscretization()
{
  using namespace chi_math::finite_element;
  chi::log.Log() << "Initializing spatial discretization.\n";
  discretization =
    chi_math::SpatialDiscretization_PWLD::New(grid, COMPUTE_CELL_MAPPINGS |
                                                    COMPUTE_UNIT_INTEGRALS);


  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log()
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi::console.GetMemoryUsageInMB() << " MB";
}