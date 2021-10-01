#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;


//###################################################################
/**Compute the steady state delayed neutron precursor concentrations.*/
void LinearBoltzmann::Solver::ComputePrecursors()
{
  auto fe =
      std::dynamic_pointer_cast<SpatialDiscretization_FE>(discretization);
  const size_t J = max_precursors_per_material;

  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  //================================================== Loop over cells
  for (const auto& cell : grid->local_cells)
  {
    const auto& fe_values = fe->GetUnitIntegrals(cell);
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double volume = transport_view.Volume();

    //==================== Obtain xs
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];

    if ((xs_id < 0) || (xs_id >= material_xs.size()))
    {
      chi_log.Log(LOG_ALLERROR)
          << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    const auto& xs = *material_xs[xs_id];

    //======================================== Loop over precursors
    for (size_t j = 0; j < xs.num_precursors; ++j)
    {
      size_t dof = cell.local_id * J + j;
      const double coeff = xs.precursor_yield[j] /
                           xs.precursor_lambda[j];

      //=================================== Loop over nodes
      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        const size_t  uk_map = transport_view.MapDOF(i, 0, 0);
        const double IntV_ShapeI = fe_values.IntV_shapeI(i);

        //============================== Loop over groups
        for (int g = 0; g < groups.size(); ++g)
          precursor_new_local[dof] += coeff * xs.nu_delayed_sigma_f[g] *
                                      phi_new_local[uk_map + g] *
                                      IntV_ShapeI / volume;
      }//for node i
    }//for precursor j

  }//for cell
}