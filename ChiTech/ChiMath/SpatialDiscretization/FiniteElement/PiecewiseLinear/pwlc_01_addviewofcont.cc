#include "pwlc.h"

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
std::shared_ptr<CellMappingFE_PWL> SpatialDiscretization_PWLC::
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
          throw std::invalid_argument("SpatialDiscretization_PWLC::MakeCellMappingFE: "
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
          throw std::invalid_argument("SpatialDiscretization_PWLC::MakeCellMappingFE: "
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
          throw std::invalid_argument("SpatialDiscretization_PWLC::MakeCellMappingFE: "
                                      "Unsupported coordinate system type encountered.");
      }
    }
    default:
      throw std::invalid_argument("SpatialDiscretization_PWLC::MakeCellPWLRawView: "
                                  "Unsupported cell type encountered.");
  }
}

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWLC::PreComputeCellSDValues()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Add cell SD-values.";

  size_t num_local_cells = ref_grid->local_cells.size();

  //================================================== Create empty view
  //                                                 for each cell
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_CELL_MAPPINGS)
    {
      chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                  << " Computing unit integrals.";
      if (!mapping_initialized)
      {
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
    chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                << " Computing unit integrals.";
      if (not integral_data_initialized)
      {
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
      chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                  << " Computing quadrature data.";
      if (not qp_data_initialized)
      {
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
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Done adding cell SD-values.";

}//AddViewOfLocalContinuum

//###################################################################
/**Returns a locally stored finite element view.*/
std::shared_ptr<CellMappingFE_PWL>
  SpatialDiscretization_PWLC::GetCellMappingFE(uint64_t cell_local_index)
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
        << "SpatialDiscretization_PWLC::MapFeView "
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
