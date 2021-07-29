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
      const auto& slab_cell = static_cast<const chi_mesh::CellSlab&>(cell);
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new SlabMappingFE_PWL(slab_cell,
                                                    ref_grid,
                                                    line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          auto cell_fe_view = new SlabMappingFE_PWL_Cylindrical(slab_cell,
                                                                ref_grid,
                                                                line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::SPHERICAL:
        {
          auto cell_fe_view = new SlabMappingFE_PWL_Spherical(slab_cell,
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
      const auto& poly_cell = static_cast<const chi_mesh::CellPolygon&>(cell);
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new PolygonMappingFE_PWL(poly_cell,
                                                       ref_grid,
                                                       tri_quad_order_arbitrary,
                                                       line_quad_order_arbitrary);

          std::shared_ptr<CellMappingFE_PWL> the_ptr(cell_fe_view);
          return the_ptr;
        }
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          auto cell_fe_view = new PolygonMappingFE_PWL_Cylindrical(poly_cell,
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
      const auto& polyh_cell = static_cast<const chi_mesh::CellPolyhedron&>(cell);
      switch (cs_type)
      {
        case chi_math::CoordinateSystemType::CARTESIAN:
        {
          auto cell_fe_view = new PolyhedronMappingFE_PWL(polyh_cell,
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
    chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                << " Computing unit integrals.";
    for (const auto& cell : ref_grid->local_cells)
      cell_mappings.push_back(MakeCellMappingFE(cell));
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                << " Computing unit integrals.";
    fe_unit_integrals.reserve(num_local_cells);
    for (size_t lc=0; lc<num_local_cells; ++lc)
    {
      UIData ui_data;

      auto cell_fe_view = GetCellMappingFE(lc);
      cell_fe_view->ComputeUnitIntegrals(ui_data);

      fe_unit_integrals.push_back(std::move(ui_data));
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
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
