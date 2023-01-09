#include "../lbts_transient_solver.h"

/**Computes the delayed neutron factor.*/
double lbs::TransientSolver::ComputeBeta()
{
  const double FPR = ComputeFissionProduction(phi_new_local);

  std::vector<double> default_zero_src(groups.size(), 0.0);

  double localDNPR = 0.0;
  for (const auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.Volume();

    auto xs = transport_view.XS();
    auto P0_src = matid_to_src_map[cell.material_id];

    const auto& S = xs.transfer_matrices;

    for (const auto& groupset : groupsets)
    {
      auto gs_i = static_cast<size_t>(groupset.groups[0].id);
      auto gs_f = static_cast<size_t>(groupset.groups.back().id);

      const auto& m_to_ell_em_map =
        groupset.quadrature->GetMomentToHarmonicsIndexMap();

      //==================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        //============================= Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          const bool fission_avail = (xs.is_fissionable and ell == 0);

          if (fission_avail)
          {
            if (options.use_precursors)
            {
              const auto& J = max_precursors_per_material;
              for (size_t j = 0; j < xs.num_precursors; ++j)
              {
                const size_t dof_map = cell.local_id * J + j;

                localDNPR += xs.chi_delayed[g][j] * xs.precursor_lambda[j] *
                             precursor_prev_local[dof_map] * cell_volume;
              }//for precursors
            }//if use precursors
          }//if fission_avail and apply_wgs_fission_src
        }//for g
      }//for m
    }//for groupset
  }//for cell

  //============================================= Allreduce global DNPR
  double DNPR = 0.0;
  MPI_Allreduce(&localDNPR, &DNPR, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return DNPR / (DNPR + FPR);
}