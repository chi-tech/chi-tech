#include "pwl.h"

#include "CHI_TECH/ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_slab.h"
#include "CHI_TECH/ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polygon.h"
#include "CHI_TECH/ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWL::PreComputeCellSDValues(
  chi_mesh::MeshContinuumPtr grid)
{
  ref_grid = grid;
  //================================================== Create empty view
  //                                                 for each cell
  if (!mapping_initialized)
  {
    cell_view_added_flags.resize(grid->local_cells.size(),false);
    mapping_initialized = true;
  }

  //================================================== Swap views for
  //                                                   specified item_id
  for (const auto& cell : grid->local_cells)
  {
    if (not cell_view_added_flags[cell.local_id])
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        auto slab_cell = (chi_mesh::CellSlab*)(&cell);
        auto cell_fe_view = new SlabPWLFEView(slab_cell,
                                              grid,
                                              line_quad_order_second);

        cell_fe_view->PreComputeValues();

        cell_fe_views.push_back(cell_fe_view);
        cell_view_added_flags[cell.local_id] = true;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto poly_cell = (chi_mesh::CellPolygon*)(&cell);
        auto cell_fe_view = new PolygonPWLFEValues(poly_cell,
                                                   grid,
                                                   tri_quad_order_second,
                                                   line_quad_order_second);

        cell_fe_view->PreComputeValues();

        cell_fe_views.push_back(cell_fe_view);
        cell_view_added_flags[cell.local_id] = true;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
        auto cell_fe_view = new PolyhedronPWLFEValues(polyh_cell,
                                                      grid,
                                                      tet_quad_order_second,
                                                      tri_quad_order_second);

        cell_fe_view->PreComputeValues();

        cell_fe_views.push_back(cell_fe_view);
        cell_view_added_flags[cell.local_id] = true;
      }
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "SpatialDiscretization_PWL::AddViewOfLocalContinuum. "
          << "Unsupported cell type encountered.";
        exit(EXIT_FAILURE);
      }
    }//if mapping not yet assigned
  }//for num cells

}//AddViewOfLocalContinuum

//###################################################################
/**Adds a PWL Finite Element for each cell of the neighboring cells.*/
void SpatialDiscretization_PWL::PreComputeNeighborCellSDValues(
  chi_mesh::MeshContinuumPtr grid)
{
  chi_log.Log(LOG_0)
    << "SpatialDiscretization_PWL::AddViewOfNeighborContinuums.";
  MPI_Barrier(MPI_COMM_WORLD);

  ref_grid = grid;

  grid->CommunicatePartitionNeighborCells(neighbor_cells);

  chi_log.Log(LOG_0)
    << "Done communicating neighbor cells.";
  MPI_Barrier(MPI_COMM_WORLD);


  //================================================== Populate cell fe views
  neighbor_cell_fe_views.reserve(neighbor_cells.size());
  for (auto cell : neighbor_cells)
  {
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      auto cell_fe_view = new SlabPWLFEView(slab_cell,
                                            grid,
                                            line_quad_order_second);

      cell_fe_view->PreComputeValues();

      neighbor_cell_fe_views.push_back(cell_fe_view);
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      auto cell_fe_view = new PolygonPWLFEValues(poly_cell,
                                                 grid,
                                                 tri_quad_order_second,
                                                 line_quad_order_second);

      cell_fe_view->PreComputeValues();

      neighbor_cell_fe_views.push_back(cell_fe_view);
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      auto cell_fe_view = new PolyhedronPWLFEValues(polyh_cell,
                                                    grid,
                                                    tet_quad_order_second,
                                                    tri_quad_order_second);

      cell_fe_view->PreComputeValues();

      neighbor_cell_fe_views.push_back(cell_fe_view);
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

}//AddViewOfNeighborContinuums


//###################################################################
/**Returns a locally stored finite element view.*/
CellPWLFEValues* SpatialDiscretization_PWL::MapFeViewL(int cell_local_index)
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

  return value;
}