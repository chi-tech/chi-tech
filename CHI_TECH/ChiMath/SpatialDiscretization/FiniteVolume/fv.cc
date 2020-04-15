#include "fv.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include "CellViews/fv_slab.h"
#include "CellViews/fv_polygon.h"
#include "CellViews/fv_polyhedron.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Only constructor for this method.*/
SpatialDiscretization_FV::SpatialDiscretization_FV(int dim)
  : SpatialDiscretization(dim)
{
  mapping_initialized = false;
}

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_FV::AddViewOfLocalContinuum(
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
      //######################################### SLAB
      if (cell.Type() == chi_mesh::CellType::SLAB)
      {
        auto view =
          new SlabFVView((chi_mesh::CellSlab*)(&cell), grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYGON
      if (cell.Type() == chi_mesh::CellType::POLYGON)
      {
        auto view =
          new PolygonFVView((chi_mesh::CellPolygon*)(&cell), grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }

      //######################################### POLYHEDRON
      if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto view =
          new PolyhedronFVView(
            (chi_mesh::CellPolyhedron*)(&cell),
            grid);

        cell_fv_views.push_back(view);
        cell_view_added_flags[cell.local_id] = true;
      }
    }//if mapping not yet assigned
  }//for num cells



}//AddViewOfLocalContinuum


//###################################################################
/**Maps the cell index to a position stored locally.*/
CellFVView* SpatialDiscretization_FV::MapFeView(int cell_local_index)
{
  CellFVView* value;
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
/**Builds finite volume based sparsity pattern.*/
void SpatialDiscretization_FV::BuildSparsityPattern(
  chi_mesh::MeshContinuum *grid,
  std::vector<int> &nodal_nnz_in_diag,
  std::vector<int> &nodal_nnz_off_diag)
{
  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = grid->local_cells.size();

  nodal_nnz_in_diag.resize(num_local_cells,0);
  nodal_nnz_off_diag.resize(num_local_cells,0);

  for (auto& cell : grid->local_cells)
  {
    int i=cell.local_id;

    nodal_nnz_in_diag[i] += 1;

    for (auto& face : cell.faces)
    {
      if (face.neighbor < 0) continue;

      if (face.IsNeighborLocal(grid))
        nodal_nnz_in_diag[i] += 1;
      else
        nodal_nnz_off_diag[i] += 1;
    }
  }
}