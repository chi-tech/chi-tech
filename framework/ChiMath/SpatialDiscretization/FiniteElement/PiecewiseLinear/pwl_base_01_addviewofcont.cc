#include "pwl_base.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void chi_math::SpatialDiscretization_PWLBase::PreComputeCellSDValues()
{
  size_t num_local_cells = ref_grid_->local_cells.size();

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags_ & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not integral_data_initialized_)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing unit integrals.";
        fe_unit_integrals_.reserve(num_local_cells);
        for (const auto& cell : ref_grid_->local_cells)
        {
          UIData ui_data;

          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.ComputeUnitIntegrals(ui_data);

          fe_unit_integrals_.push_back(std::move(ui_data));
        }

        integral_data_initialized_ = true;
      }
    }//if compute unit intgrls
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags_ & SetupFlags::COMPUTE_QP_DATA)
    {
      if (not qp_data_initialized_)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing quadrature data.";
        fe_vol_qp_data_.reserve(num_local_cells);
        fe_srf_qp_data_.reserve(num_local_cells);
        for (const auto& cell : ref_grid_->local_cells)
        {
          fe_vol_qp_data_.emplace_back();
          fe_srf_qp_data_.emplace_back();

          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeAllQuadraturePointData(fe_vol_qp_data_.back(),
                                                        fe_srf_qp_data_.back());
        }

        qp_data_initialized_ = true;
      }
    }//if init qp data
  }

}//AddViewOfLocalContinuum


//###################################################################
/**Adds a PWL Finite Element for each cell of the neighboring cells.*/
void chi_math::SpatialDiscretization_PWLBase::PreComputeNeighborCellSDValues()
{
  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags_ & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not nb_integral_data_initialized_)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing neighbor unit integrals.";
        const auto ghost_ids = ref_grid_->cells.GetGhostGlobalIDs();
        for (uint64_t cell_gid : ghost_ids)
        {
          const auto& cell = ref_grid_->cells[cell_gid];
          const auto& cell_mapping = GetCellMapping(cell);

          UIData ui_data;
          cell_mapping.ComputeUnitIntegrals(ui_data);

          using namespace std;
          nb_fe_unit_integrals_.insert(make_pair(cell_gid, move(ui_data)));
        }

        nb_integral_data_initialized_ = true;
      }
    }//if compute unit intgrls
  }

  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags_ & SetupFlags::COMPUTE_QP_DATA)
    {
      if (not nb_qp_data_initialized_)
      {
        chi::log.Log() << chi::program_timer.GetTimeString()
                       << " Computing neighbor quadrature data.";
        const auto ghost_ids = ref_grid_->cells.GetGhostGlobalIDs();
        for (uint64_t cell_gid : ghost_ids)
        {
          const auto& cell = ref_grid_->cells[cell_gid];
          const auto& cell_mapping = GetCellMapping(cell);

          QPDataVol qp_data_vol;
          std::vector<QPDataFace> qp_data_srf;
          cell_mapping.InitializeAllQuadraturePointData(qp_data_vol,
                                                        qp_data_srf);

          using namespace std;
          nb_fe_vol_qp_data_.insert(make_pair(cell_gid, move(qp_data_vol)));
          nb_fe_srf_qp_data_.insert(make_pair(cell_gid, move(qp_data_srf)));
        }

        nb_qp_data_initialized_ = true;
      }
    }//if init qp data
  }

}//AddViewOfNeighborContinuums