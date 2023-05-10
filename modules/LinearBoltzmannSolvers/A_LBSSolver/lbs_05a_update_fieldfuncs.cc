#include "lbs_solver.h"

#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_log.h"

// ###################################################################
/**Copy relevant section of phi_old to the field functions.*/
void lbs::LBSSolver::UpdateFieldFunctions()
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  //======================================== Update flux moments
  for (size_t g = 0; g < groups_.size(); ++g)
  {
    for (size_t m = 0; m < num_moments_; ++m)
    {
      std::vector<double> data_vector_local(local_node_count_, 0.0);

      for (const auto& cell : grid_ptr_->local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.NumNodes();

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          const int64_t imapB = sdm.MapDOFLocal(cell, i);

          data_vector_local[imapB] = phi_old_local_[imapA];
        } // for node
      }   // for cell

      ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
                        "Update error for phi based field functions");

      const size_t ff_index = phi_field_functions_local_map_.at({g, m});

      auto& ff_ptr = field_functions_.at(ff_index);
      ff_ptr->UpdateFieldVector(data_vector_local);
    } // for m
  }   // for g

  //======================================== Update power generation
  //                                         if enabled
  if (options_.power_field_function_on)
  {
    std::vector<double> data_vector_local(local_node_count_, 0.0);

    double local_total_power = 0.0;
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      const auto& Vi = unit_cell_matrices_[cell.local_id_].Vi_vectors;

      const auto& xs = matid_to_xs_map_.at(cell.material_id_);

      if (not xs->IsFissionable()) continue;

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i);
        const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

        double nodal_power = 0.0;
        for (size_t g = 0; g < groups_.size(); ++g)
        {
          const double sigma_fg = xs->SigmaFission()[g];
          // const double kappa_g = xs->Kappa()[g];
          const double kappa_g = options_.power_default_kappa;

          nodal_power += kappa_g * sigma_fg * phi_old_local_[imapB + g];
        } // for g

        data_vector_local[imapA] = nodal_power;
        local_total_power += nodal_power * Vi[i];
      } // for node
    }   // for cell

    if (options_.power_normalization > 0.0)
    {
      double globl_total_power;
      MPI_Allreduce(&local_total_power, // sendbuf
                    &globl_total_power, // recvbuf
                    1,
                    MPI_DOUBLE,      // count + datatype
                    MPI_SUM,         // operation
                    MPI_COMM_WORLD); // communicator

      chi_math::Scale(data_vector_local,
                      options_.power_normalization / globl_total_power);
    }

    const size_t ff_index = power_gen_fieldfunc_local_handle_;

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);

  } // if power enabled
}
