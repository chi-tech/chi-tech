#include "pwlc.h"

#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void chi_math::SpatialDiscretization_PWLC::PreComputeCellSDValues()
{
  size_t num_local_cells = ref_grid->local_cells.size();

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not integral_data_initialized)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing unit integrals.";
        fe_unit_integrals.reserve(num_local_cells);
        for (const auto& cell : ref_grid->local_cells)
        {
          UIData ui_data;

          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.ComputeUnitIntegrals(ui_data);

          fe_unit_integrals.push_back(std::move(ui_data));
        }

        integral_data_initialized = true;
      }
    }//if compute unit intgrls
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_QP_DATA)
    {
      if (not qp_data_initialized)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing quadrature data.";
        fe_vol_qp_data.reserve(num_local_cells);
        fe_srf_qp_data.reserve(num_local_cells);
        for (const auto& cell : ref_grid->local_cells)
        {
          fe_vol_qp_data.emplace_back();
          fe_srf_qp_data.emplace_back();

          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeAllQuadraturePointData(fe_vol_qp_data.back(),
                                                        fe_srf_qp_data.back());
        }

        qp_data_initialized = true;
      }
    }//if init qp data
  }

}//AddViewOfLocalContinuum
