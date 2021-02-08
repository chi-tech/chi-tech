#include "fv.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include "CellViews/fv_slab.h"
#include "CellViews/fv_polygon.h"
#include "CellViews/fv_polyhedron.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_FV::PreComputeCellSDValues(
  chi_mesh::MeshContinuumPtr grid)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "SpatialDiscretization_FV - Adding view of local continuum.";
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
      //######################################### SLAB
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        auto view =
          new SlabFVValues((chi_mesh::CellSlab*)(&cell), grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYGON
      if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto view =
          new PolygonFVValues((chi_mesh::CellPolygon*)(&cell), grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYHEDRON
      if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto view =
          new PolyhedronFVValues(
            (chi_mesh::CellPolyhedron*)(&cell),
            grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }
    }//if mapping not yet assigned
  }//for num cells
}//AddViewOfLocalContinuum

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_FV::
  AddViewOfNeighborContinuums(chi_mesh::MeshContinuumPtr grid)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "SpatialDiscretization_FV - Adding view of neighbor continuums.";
  ref_grid = grid;

  grid->CommunicatePartitionNeighborCells(neighbor_cells);

  chi_log.Log(LOG_0VERBOSE_1)
    << "Number neighbor cells: " << neighbor_cells.size();

  //================================================== Reorder according to
  //                                                   ghost indices
  std::vector<chi_mesh::Cell*> temp(grid->cells.GetNumGhosts(), nullptr);

  auto ghost_ids = grid->cells.GetGhostGlobalIDs();
  chi_log.Log(LOG_0VERBOSE_1)
    << "Number of ghost ids: " << ghost_ids.size();
  for (int ghost_id : ghost_ids)
  {
    int ghost_local_index = grid->cells.GetGhostLocalID(ghost_id);

    for (auto cell : neighbor_cells)
      if (cell->global_id == ghost_id)
        temp[ghost_local_index] = cell;
  }
  neighbor_cells.clear();
  for (auto tcell : temp)
    if (tcell != nullptr)
      neighbor_cells.push_back(tcell);


  chi_log.Log(LOG_0VERBOSE_1)
    << "Adding neighbor views";
  //================================================== Populate cell fe views
  neighbor_cell_fv_views.reserve(neighbor_cells.size());
  for (auto cell : neighbor_cells)
  {
    //######################################### SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto view =
        new SlabFVValues((chi_mesh::CellSlab*)cell, grid);

      neighbor_cell_fv_views.push_back(view);
    }

    //######################################### POLYGON
    if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto view =
        new PolygonFVValues((chi_mesh::CellPolygon*)cell, grid);

      neighbor_cell_fv_views.push_back(view);
    }

    //######################################### POLYHEDRON
    if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto view =
        new PolyhedronFVValues(
          (chi_mesh::CellPolyhedron*)cell,
          grid);

      neighbor_cell_fv_views.push_back(view);
    }
  }//for num cells


  chi_log.Log(LOG_ALLVERBOSE_1)
    << "Number of neighbor cells added: "
    << neighbor_cell_fv_views.size();
}//AddViewOfNeighborContinuums

//###################################################################
/**Maps the cell index to a position stored locally.*/
CellFVValues* SpatialDiscretization_FV::MapFeView(int cell_local_index)
{
  CellFVValues* value;
  try { value = cell_fv_views.at(cell_local_index); }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapFeView "
         "Failure to map Finite Volume View. The view is either not"
         "available or the supplied local index is invalid.";
    exit(EXIT_FAILURE);
  }

  return value;
}

//###################################################################
/**Maps the cell index to a position stored locally.*/
CellFVValues* SpatialDiscretization_FV::MapNeighborFeView(int cell_global_index)
{
  auto& cell = ref_grid->cells[cell_global_index];

  if (cell.partition_id == chi_mpi.location_id)
    return MapFeView(cell.local_id);
  else
  {
    int index=0;
    for (auto ncell : neighbor_cells)
    {
      if (ncell->global_id == cell_global_index)
        break;
      ++index;
    }
    return neighbor_cell_fv_views[index];
  }
}