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
    cell_fe_views_mapping.resize(grid->cells.size(), -1);
    mapping_initialized = true;
  }


  //================================================== Swap views for
  //                                                   specified item_id
  for (auto& cell_index : grid->local_cell_glob_indices)
  {
    chi_mesh::Cell* cell = grid->cells[cell_index];

    if (cell_fe_views_mapping[cell_index]<0)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell->Type() == chi_mesh::CellType::SLAB)
      {
        auto slab_cell = dynamic_cast<chi_mesh::CellSlab*>(cell);
        auto cell_fe_view = new SlabFEView(slab_cell, grid);

        //cell_fe_view->PreCompute();

        this->cell_fe_views.push_back(cell_fe_view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell->Type() == chi_mesh::CellType::POLYGON)
      {
        auto poly_cell = dynamic_cast<chi_mesh::CellPolygon*>(cell);
        auto cell_fe_view = new PolygonFEView(poly_cell, grid, this);

        cell_fe_view->PreCompute();

        this->cell_fe_views.push_back(cell_fe_view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto polyh_cell = dynamic_cast<chi_mesh::CellPolyhedron*>(cell);
        auto cell_fe_view = new PolyhedronFEView(polyh_cell, grid, this);

        cell_fe_view->PreCompute();
        cell_fe_view->CleanUp();
        this->cell_fe_views.push_back(cell_fe_view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
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


  chi_log.Log(LOG_ALL)
    << "Number of neighbor cells added: "
    << neighbor_cell_fe_views.size();

}//AddViewOfNeighborContinuums






//###################################################################
/**Maps the cell index to a position stored locally.*/
CellFEView* SpatialDiscretization_PWL::MapFeView(int cell_glob_index)
{
  CellFEView* value = cell_fe_views.at(cell_fe_views_mapping[cell_glob_index]);
  return value;
}