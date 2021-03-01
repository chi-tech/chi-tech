#include "k_eigenvalue_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

using namespace LinearBoltzmann;

//###################################################################
/**Compute the volume integral of vector phi.*/
std::vector<double> 
KEigenvalue::Solver::IntegrateVolume(std::vector<double> phi)
{
  std::vector<double> Phi(num_groups,0.0);

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  // ----- Loop over cells
  for (auto& cell : grid->local_cells)
  {
    // ----- Cell information
    auto cell_fe_view = pwl->GetUnitIntegrals(cell);
    auto transport_view = cell_transport_views[cell.local_id];

    // ----- Loop over cell dofs
    for (int i=0; i<cell_fe_view.NumNodes(); i++)
    {
      int ir = transport_view.MapDOF(i,0,0);
      double* phip = &phi[ir];

      double intV_shapeI = cell_fe_view.IntV_shapeI(i);

      // ----- Loop over groups and integrate
      for (int g=first_grp; g<=last_grp; g++)
        Phi[g] += phip[g] * intV_shapeI;      
    }//for i
  }//for c
  return Phi;
}