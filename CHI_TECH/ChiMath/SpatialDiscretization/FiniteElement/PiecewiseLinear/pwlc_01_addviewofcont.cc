#include "pwlc.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWLC::PreComputeCellSDValues()
{
  //================================================== Create empty view
  //                                                 for each cell
  if (!mapping_initialized)
  {
    for (const auto& cell : ref_grid->local_cells)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        auto slab_cell = (chi_mesh::CellSlab*)(&cell);
        auto cell_fe_view = new SlabPWLFEView(slab_cell,
                                              ref_grid,
                                              line_quad_order_second,
                                              line_quad_order_arbitrary);

        cell_fe_views.push_back(cell_fe_view);
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto poly_cell = (chi_mesh::CellPolygon*)(&cell);
        auto cell_fe_view = new PolygonPWLFEValues(poly_cell,
                                                   ref_grid,
                                                   tri_quad_order_second,
                                                   line_quad_order_second,
                                                   tri_quad_order_arbitrary,
                                                   line_quad_order_arbitrary);

        cell_fe_views.push_back(cell_fe_view);
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
        auto cell_fe_view = new PolyhedronPWLFEValues(polyh_cell,
                                                      ref_grid,
                                                      tet_quad_order_second,
                                                      tri_quad_order_second,
                                                      tet_quad_order_arbitrary,
                                                      tri_quad_order_arbitrary);

        cell_fe_views.push_back(cell_fe_view);
      }
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "SpatialDiscretization_PWL::AddViewOfLocalContinuum. "
          << "Unsupported cell type encountered.";
        exit(EXIT_FAILURE);
      }
    }//for num cells

    mapping_initialized = true;
  }

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not integral_data_initialized)
      {
        fe_unit_integrals.reserve(cell_fe_views.size());
        for (auto& cell_fe_view : cell_fe_views)
        {
          UIData ui_data;
          cell_fe_view->ComputeUnitIntegrals(ui_data);

          fe_unit_integrals.push_back(std::move(ui_data));
        }

        integral_data_initialized = true;
      }
    }//if compute unit intgrls
  }


  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::INIT_QP_DATA)
    {
      if (not qp_data_initialized)
      {
        fe_vol_qp_data.reserve(cell_fe_views.size());
        fe_srf_qp_data.reserve(cell_fe_views.size());
        for (auto& cell_fe_view : cell_fe_views)
        {
          fe_vol_qp_data.emplace_back();
          fe_srf_qp_data.emplace_back();
          cell_fe_view->InitializeQuadraturePointData(fe_vol_qp_data.back(),
                                                      fe_srf_qp_data.back());
        }

        qp_data_initialized = true;
      }
    }//if init qp data
  }


}//AddViewOfLocalContinuum

//###################################################################
/**Returns a locally stored finite element view.*/
CellPWLFEValues& SpatialDiscretization_PWLC::GetCellFEView(int cell_local_index)
{
  CellPWLFEValues* value;
  try { value = cell_fe_views.at(cell_local_index); }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_PWL::MapFeView "
         "Failure to map Finite Element View. The view is either not"
         "available or the supplied local index is invalid.";
    exit(EXIT_FAILURE);
  }

  return *value;
}