#include "pwlc.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Makes a shared_ptr CellPWLView for a cell based on its type.*/
std::shared_ptr<CellPWLFEValues> SpatialDiscretization_PWLC::
  MakeCellPWLView(const chi_mesh::Cell& cell) const
{
  switch (cell.Type())
  {
    case chi_mesh::CellType::SLAB:
    {
      const auto& slab_cell = (const chi_mesh::CellSlab&)(cell);
      auto cell_fe_view = new SlabPWLFEView(slab_cell,
                                            ref_grid,
                                            line_quad_order_second,
                                            line_quad_order_arbitrary);

      std::shared_ptr<CellPWLFEValues> the_ptr(cell_fe_view);
      return the_ptr;
    }
    case chi_mesh::CellType::POLYGON:
    {
      const auto& poly_cell = (const chi_mesh::CellPolygon&)(cell);
      auto cell_fe_view = new PolygonPWLFEValues(poly_cell,
                                                 ref_grid,
                                                 tri_quad_order_second,
                                                 line_quad_order_second,
                                                 tri_quad_order_arbitrary,
                                                 line_quad_order_arbitrary);

      std::shared_ptr<CellPWLFEValues> the_ptr(cell_fe_view);
      return the_ptr;
    }
    case chi_mesh::CellType::POLYHEDRON:
    {
      const auto& polyh_cell = (const chi_mesh::CellPolyhedron&)(cell);
      auto cell_fe_view = new PolyhedronPWLFEValues(polyh_cell,
                                                    ref_grid,
                                                    tet_quad_order_second,
                                                    tri_quad_order_second,
                                                    tet_quad_order_arbitrary,
                                                    tri_quad_order_arbitrary);

      std::shared_ptr<CellPWLFEValues> the_ptr(cell_fe_view);
      return the_ptr;
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
    if (setup_flags & SetupFlags::COMPUTE_CELL_VIEWS)
    {
      chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                                  << " Computing unit integrals.";
      if (!mapping_initialized)
      {
        for (const auto& cell : ref_grid->local_cells)
          cell_fe_views.push_back(MakeCellPWLView(cell));

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

          auto cell_fe_view = GetCellPWLView(lc);
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
    if (setup_flags & SetupFlags::INIT_QP_DATA)
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

          auto cell_fe_view = GetCellPWLView(lc);
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
std::shared_ptr<CellPWLFEValues>
  SpatialDiscretization_PWLC::GetCellPWLView(int cell_local_index)
{

  if (mapping_initialized)
  {
    try
    {
      return cell_fe_views.at(cell_local_index);
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
    return MakeCellPWLView(ref_grid->local_cells[cell_local_index]);
  }

}