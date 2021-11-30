#include "pwl.h"

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polyhedron.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_cylindrical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_spherical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon_cylindrical.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Makes a shared_ptr CellPWLView for a cell based on its type.*/
std::shared_ptr<CellMappingFE_PWL> SpatialDiscretization_PWLD::
  MakeCellMappingFE(const chi_mesh::Cell& cell) const
{
  switch (cell.Type())
  {
    case chi_mesh::CellType::SLAB:
    {
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new SlabMappingFE_PWL(cell,
                                                    ref_grid,
                                                    line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          auto cell_fe_view = new SlabMappingFE_PWL_Cylindrical(cell,
                                                                ref_grid,
                                                                line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::SPHERICAL:
        {
          auto cell_fe_view = new SlabMappingFE_PWL_Spherical(cell,
                                                              ref_grid,
                                                              line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        default:
          throw std::invalid_argument("SpatialDiscretization_PWL::MakeCellMappingFE: "
                                      "Unsupported coordinate system type encountered.");
      }
    }
    case chi_mesh::CellType::POLYGON:
    {
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new PolygonMappingFE_PWL(cell,
                                                       ref_grid,
                                                       tri_quad_order_arbitrary,
                                                       line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          auto cell_fe_view = new PolygonMappingFE_PWL_Cylindrical(cell,
                                                                   ref_grid,
                                                                   tri_quad_order_arbitrary,
                                                                   line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        default:
          throw std::invalid_argument("SpatialDiscretization_PWL::MakeCellMappingFE: "
                                      "Unsupported coordinate system type encountered.");
      }
    }
    case chi_mesh::CellType::POLYHEDRON:
    {
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new PolyhedronMappingFE_PWL(cell,
                                                          ref_grid,
                                                          tet_quad_order_arbitrary,
                                                          tri_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        default:
          throw std::invalid_argument("SpatialDiscretization_PWL::MakeCellMappingFE: "
                                      "Unsupported coordinate system type encountered.");
      }
    }
    default:
      throw std::invalid_argument("SpatialDiscretization_PWL::MakeCellPWLRawView: "
                                  "Unsupported cell type encountered.");
  }
}

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWLD::PreComputeCellSDValues()
{
  size_t num_local_cells = ref_grid->local_cells.size();

  //================================================== Create empty view
  //                                                 for each cell
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_CELL_MAPPINGS)
    {
      if (!mapping_initialized)
      {
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing cell views";
        for (const auto& cell : ref_grid->local_cells)
          cell_mappings.push_back(MakeCellMappingFE(cell));

        mapping_initialized = true;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not integral_data_initialized)
      {
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing unit integrals.";
        fe_unit_integrals.reserve(num_local_cells);
        for (size_t lc=0; lc<num_local_cells; ++lc)
        {
          UIData ui_data;

          auto cell_fe_view = GetCellMappingFE(lc);
          cell_fe_view->ComputeUnitIntegrals(ui_data);

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
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing quadrature data.";
        fe_vol_qp_data.reserve(num_local_cells);
        fe_srf_qp_data.reserve(num_local_cells);
        for (size_t lc=0; lc<num_local_cells; ++lc)
        {
          fe_vol_qp_data.emplace_back();
          fe_srf_qp_data.emplace_back();

          auto cell_fe_view = GetCellMappingFE(lc);
          cell_fe_view->InitializeAllQuadraturePointData(fe_vol_qp_data.back(),
                                                         fe_srf_qp_data.back());
        }

        qp_data_initialized = true;
      }
    }//if init qp data
  }
}//AddViewOfLocalContinuum

//###################################################################
/**Adds a PWL Finite Element for each cell of the neighboring cells.*/
void SpatialDiscretization_PWLD::PreComputeNeighborCellSDValues()
{
  //================================================== Populate cell fe views
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_CELL_MAPPINGS)
    {
      if (not nb_mapping_initialized)
      {
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing neighbor cell views.";
        for (auto& cell_map : neighbor_cells)
        {
          const auto& cell = *cell_map.second;

          neighbor_cell_fe_views.insert(
            std::make_pair(cell.global_id, MakeCellMappingFE(cell)));
        }//for num cells

        nb_mapping_initialized = true;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not nb_integral_data_initialized)
      {
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing neighbor unit integrals.";
        for (auto& nb_cell : neighbor_cells)
        {
          uint64_t cell_global_id = nb_cell.first;
          auto cell_fe_view = GetNeighborCellMappingFE(cell_global_id);

          UIData ui_data;
          cell_fe_view->ComputeUnitIntegrals(ui_data);

          nb_fe_unit_integrals.insert(std::make_pair(cell_global_id,std::move(ui_data)));
        }

        nb_integral_data_initialized = true;
      }
    }//if compute unit intgrls
  }


  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_QP_DATA)
    {
      if (not nb_qp_data_initialized)
      {
        chi_log.Log() << chi_program_timer.GetTimeString()
                      << " Computing neighbor quadrature data.";
        for (auto& nb_cell : neighbor_cells)
        {
          uint64_t cell_global_id = nb_cell.first;
          auto cell_fe_view = GetNeighborCellMappingFE(cell_global_id);

          QPDataVol qp_data_vol;
          std::vector<QPDataFace> qp_data_srf;
          cell_fe_view->InitializeAllQuadraturePointData(qp_data_vol,
                                                         qp_data_srf);

          nb_fe_vol_qp_data.insert(std::make_pair(cell_global_id,std::move(qp_data_vol)));
          nb_fe_srf_qp_data.insert(std::make_pair(cell_global_id,std::move(qp_data_srf)));
        }

        nb_qp_data_initialized = true;
      }
    }//if init qp data
  }

}//AddViewOfNeighborContinuums


//###################################################################
/**Returns a locally stored finite element view.*/
std::shared_ptr<CellMappingFE_PWL>
  SpatialDiscretization_PWLD::GetCellMappingFE(uint64_t cell_local_index)
{
  if (mapping_initialized)
  {
    try
    {
      return cell_mappings.at(cell_local_index);
    }
    catch (const std::out_of_range& o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "SpatialDiscretization_PWL::MapFeView "
           "Failure to map Finite Element View. The view is either not"
           "available or the supplied local index is invalid.";
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    return MakeCellMappingFE(ref_grid->local_cells[cell_local_index]);
  }
}

//###################################################################
/**Maps a neigboring cell from a global cell index. The spatial discretizations
 * maintains a non-ghost version of all neighboring cells.*/
chi_mesh::Cell& SpatialDiscretization_PWLD::
  GetNeighborCell(uint64_t cell_glob_index)
{
  //=================================== First check locally
  if (ref_grid->IsCellLocal(cell_glob_index))
    return ref_grid->cells[cell_glob_index];

  //=================================== Now check neighbor cells
  auto neighbor_location = neighbor_cells.find(cell_glob_index);

  if (neighbor_location != neighbor_cells.end())
    return *neighbor_cells.at(cell_glob_index);
  else
    throw std::logic_error(std::string(__FUNCTION__) +
                           " Mapping of neighbor cell failed.");
}

//###################################################################
/**Maps a neigboring cell's fe view from a global cell index.*/
std::shared_ptr<CellMappingFE_PWL> SpatialDiscretization_PWLD::
  GetNeighborCellMappingFE(uint64_t cell_glob_index)
{
  //=================================== First check locally
  if (ref_grid->IsCellLocal(cell_glob_index))
  {
    auto& neighbor_cell = ref_grid->cells[cell_glob_index];
    return GetCellMappingFE(neighbor_cell.local_id);
  }

  //=================================== Now check neighbor cells
  if (nb_mapping_initialized)
  {
    auto neighbor_location = neighbor_cell_fe_views.find(cell_glob_index);

    if (neighbor_location != neighbor_cell_fe_views.end())
      return neighbor_cell_fe_views.at(cell_glob_index);
    else
      throw std::logic_error(std::string(__FUNCTION__) +
                             " Mapping of neighbor cell failed.");
  }
  else
  {
    return MakeCellMappingFE(GetNeighborCell(cell_glob_index));
  }

}
