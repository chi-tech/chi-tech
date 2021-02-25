#include "pwl.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWL::PreComputeCellSDValues()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Add cell SD-values.";

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
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0VERBOSE_1) << chi_program_timer.GetTimeString()
                              << " Done adding cell SD-values.";

}//AddViewOfLocalContinuum

//###################################################################
/**Adds a PWL Finite Element for each cell of the neighboring cells.*/
void SpatialDiscretization_PWL::PreComputeNeighborCellSDValues()
{
  chi_log.Log(LOG_0)
    << "SpatialDiscretization_PWL::AddViewOfNeighborContinuums.";
  MPI_Barrier(MPI_COMM_WORLD);

  ref_grid->CommunicatePartitionNeighborCells(neighbor_cells);

  chi_log.Log(LOG_0)
    << "Done communicating neighbor cells.";
  MPI_Barrier(MPI_COMM_WORLD);


  //================================================== Populate cell fe views
  for (auto& cell_map : neighbor_cells)
  {
    auto& cell = cell_map.second;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      auto cell_fe_view = new SlabPWLFEView(slab_cell,
                                            ref_grid,
                                            line_quad_order_second,
                                            line_quad_order_arbitrary);

      neighbor_cell_fe_views.insert(std::pair<uint64_t, CellPWLFEValues*>(
        cell->global_id,cell_fe_view));
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      auto cell_fe_view = new PolygonPWLFEValues(poly_cell,
                                                 ref_grid,
                                                 tri_quad_order_second,
                                                 line_quad_order_second,
                                                 tri_quad_order_arbitrary,
                                                 line_quad_order_arbitrary);

      neighbor_cell_fe_views.insert(std::pair<uint64_t, CellPWLFEValues*>(
        cell->global_id,cell_fe_view));
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      auto cell_fe_view = new PolyhedronPWLFEValues(polyh_cell,
                                                    ref_grid,
                                                    tet_quad_order_second,
                                                    tri_quad_order_second,
                                                    tet_quad_order_arbitrary,
                                                    tri_quad_order_arbitrary);

      neighbor_cell_fe_views.insert(std::pair<uint64_t, CellPWLFEValues*>(
        cell->global_id,cell_fe_view));
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "SpatialDiscretization_PWL::AddViewOfNeighborContinuums. "
        << "Unsupported cell type encountered.";
      exit(EXIT_FAILURE);
    }
  }//for num cells


  chi_log.Log(LOG_ALLVERBOSE_1)
    << "Number of neighbor cells added: "
    << neighbor_cell_fe_views.size();

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
      if (not nb_integral_data_initialized)
      {
        for (auto& cell_fe_view_pair : neighbor_cell_fe_views)
        {
          uint64_t cell_global_id = cell_fe_view_pair.first;
          auto cell_fe_view = cell_fe_view_pair.second;

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
    if (setup_flags & SetupFlags::INIT_QP_DATA)
    {
      if (not nb_qp_data_initialized)
      {
        for (auto& cell_fe_view_pair : neighbor_cell_fe_views)
        {
          uint64_t cell_global_id = cell_fe_view_pair.first;
          auto cell_fe_view = cell_fe_view_pair.second;

          QPDataVol qp_data_vol;
          std::vector<QPDataFace> qp_data_srf;
          cell_fe_view->InitializeQuadraturePointData(qp_data_vol,
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
CellPWLFEValues& SpatialDiscretization_PWL::GetCellFEView(int cell_local_index)
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

//###################################################################
/**Maps a neigboring cell from a global cell index.*/
chi_mesh::Cell& SpatialDiscretization_PWL::
  MapNeighborCell(int cell_glob_index)
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
CellPWLFEValues& SpatialDiscretization_PWL::
  GetNeighborCellFEView(int cell_glob_index)
{
  //=================================== First check locally
  if (ref_grid->IsCellLocal(cell_glob_index))
  {
    auto& neighbor_cell = ref_grid->cells[cell_glob_index];
    return GetCellFEView(neighbor_cell.local_id);
  }

  //=================================== Now check neighbor cells
  auto neighbor_location = neighbor_cell_fe_views.find(cell_glob_index);

  if (neighbor_location != neighbor_cell_fe_views.end())
    return *neighbor_cell_fe_views.at(cell_glob_index);
  else
    throw std::logic_error(std::string(__FUNCTION__) +
                           " Mapping of neighbor cell failed.");
}