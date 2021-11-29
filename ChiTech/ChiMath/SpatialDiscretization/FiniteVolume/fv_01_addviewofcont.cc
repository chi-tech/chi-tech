#include "fv.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include "CellViews/fv_slab.h"
#include "CellViews/fv_polygon.h"
#include "CellViews/fv_polyhedron.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_FV::PreComputeCellSDValues()
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "SpatialDiscretization_FV - Adding view of local continuum.";

  //================================================== Create empty view
  //                                                 for each cell
  if (!mapping_initialized)
  {
    cell_view_added_flags.resize(ref_grid->local_cells.size(),false);
    mapping_initialized = true;
  }


  //================================================== Swap views for
  //                                                   specified item_id
  for (const auto& cell : ref_grid->local_cells)
  {
    if (not cell_view_added_flags[cell.local_id])
    {
      //######################################### SLAB
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        auto view = new SlabFVValues(cell, *ref_grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYGON
      if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto view = new PolygonFVValues(cell, *ref_grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYHEDRON
      if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto view = new PolyhedronFVValues(cell, *ref_grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }
    }//if mapping not yet assigned
  }//for num cells
}//AddViewOfLocalContinuum

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_FV::
  PreComputeNeighborCellSDValues()
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "SpatialDiscretization_FV - Adding view of neighbor continuums.";

  ref_grid->CommunicatePartitionNeighborCells(neighbor_cells);

  chi_log.Log(LOG_0VERBOSE_1)
    << "Number neighbor cells: " << neighbor_cells.size();

  //================================================== Reorder according to
  //                                                   ghost indices
  std::vector<chi_mesh::Cell*> temp(ref_grid->cells.GetNumGhosts(), nullptr);

  auto ghost_ids = ref_grid->cells.GetGhostGlobalIDs();
  chi_log.Log(LOG_0VERBOSE_1)
    << "Number of ghost ids: " << ghost_ids.size();
  for (uint64_t ghost_id : ghost_ids)
  {
    uint64_t ghost_local_index = ref_grid->cells.GetGhostLocalID(ghost_id);

    for (auto& cell_pair : neighbor_cells)
      if (cell_pair.first == ghost_id)
        temp[ghost_local_index] = cell_pair.second;
  }
  neighbor_cells.clear();
  for (auto tcell : temp)
    if (tcell != nullptr)
      neighbor_cells.insert(std::pair<uint64_t, chi_mesh::Cell*>(
        tcell->global_id, tcell));


  chi_log.Log(LOG_0VERBOSE_1)
    << "Adding neighbor views";
  //================================================== Populate cell fe views
  for (auto& cell_pair : neighbor_cells)
  {
    const auto& cell = *cell_pair.second;
    //######################################### SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto view = new SlabFVValues(cell, *ref_grid);

      neighbor_cell_fv_views.insert(std::pair<uint64_t, CellFVValues*>(
        cell.global_id,view));
    }

    //######################################### POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto view = new PolygonFVValues(cell, *ref_grid);

      neighbor_cell_fv_views.insert(std::pair<uint64_t, CellFVValues*>(
        cell.global_id,view));
    }

    //######################################### POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto view = new PolyhedronFVValues(cell, *ref_grid);

      neighbor_cell_fv_views.insert(std::pair<uint64_t, CellFVValues*>(
        cell.global_id,view));
    }
  }//for num cells


  chi_log.Log(LOG_ALLVERBOSE_1)
    << "Number of neighbor cells added: "
    << neighbor_cell_fv_views.size();
}//AddViewOfNeighborContinuums

//###################################################################
/**Maps the cell index to a position stored locally.*/
CellFVValues* SpatialDiscretization_FV::MapFeView(uint64_t cell_local_index)
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
CellFVValues* SpatialDiscretization_FV::MapNeighborFeView(uint64_t cell_global_index)
{
  //=================================== First check locally
  if (ref_grid->IsCellLocal(cell_global_index))
  {
    auto& neighbor_cell = ref_grid->cells[cell_global_index];
    return MapFeView(neighbor_cell.local_id);
  }

  //=================================== Now check neighbor cells
  auto neighbor_location = neighbor_cell_fv_views.find(cell_global_index);

  if (neighbor_location != neighbor_cell_fv_views.end())
    return neighbor_cell_fv_views.at(cell_global_index);
  else
    throw std::logic_error(std::string(__FUNCTION__) +
                           " Mapping of neighbor cell failed.");


//  auto& cell = ref_grid->cells[cell_global_index];
//
//  if (cell.partition_id == chi_mpi.location_id)
//    return MapFeView(cell.local_id);
//  else
//  {
//    int index=0;
//    for (auto ncell : neighbor_cells)
//    {
//      if (ncell->global_id == cell_global_index)
//        break;
//      ++index;
//    }
//    return neighbor_cell_fv_views[index];
//  }
}