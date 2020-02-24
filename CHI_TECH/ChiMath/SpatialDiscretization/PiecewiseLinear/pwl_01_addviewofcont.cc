#include "pwl.h"

#include "CellViews/pwl_slab.h"
#include "CellViews/pwl_polygon.h"
#include "CellViews/pwl_polyhedron.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWL::AddViewOfLocalContinuum(
  chi_mesh::MeshContinuum* grid)
{
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
        auto cell_fe_view = new SlabFEView(slab_cell, grid);

        //cell_fe_view->PreCompute();
        //cell_fe_view->CleanUp();
        cell_fe_views.push_back(cell_fe_view);
        cell_view_added_flags[cell.local_id] = true;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto poly_cell = (chi_mesh::CellPolygon*)(&cell);
        auto cell_fe_view = new PolygonFEView(poly_cell, grid, this);

        cell_fe_view->PreCompute();

        cell_fe_views.push_back(cell_fe_view);
        cell_view_added_flags[cell.local_id] = true;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
        auto cell_fe_view = new PolyhedronFEView(polyh_cell, grid, this);

        cell_fe_view->PreCompute();
        cell_fe_view->CleanUp();
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
void SpatialDiscretization_PWL::AddViewOfNeighborContinuums(
  chi_mesh::MeshContinuum* grid)
{
  grid->CommunicatePartitionNeighborCells(neighbor_cells);


  //================================================== Populate cell fe views
  neighbor_cell_fe_views.reserve(neighbor_cells.size());
  for (auto cell : neighbor_cells)
  {
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      auto cell_fe_view = new SlabFEView(slab_cell, grid);

      //cell_fe_view->PreCompute();

      neighbor_cell_fe_views.push_back(cell_fe_view);
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      auto cell_fe_view = new PolygonFEView(poly_cell, grid, this);

      cell_fe_view->PreCompute();

      neighbor_cell_fe_views.push_back(cell_fe_view);
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      auto cell_fe_view = new PolyhedronFEView(polyh_cell, grid, this);

      cell_fe_view->PreCompute();
      cell_fe_view->CleanUp();
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
CellFEView* SpatialDiscretization_PWL::MapFeViewL(int cell_local_index)
{
  CellFEView* value;
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