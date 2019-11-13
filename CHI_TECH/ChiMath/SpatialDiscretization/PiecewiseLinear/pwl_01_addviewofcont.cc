#include "pwl.h"

#include "CellViews/pwl_slab.h"
#include "CellViews/pwl_polygon.h"
#include "CellViews/pwl_polyhedron.h"

#include <chi_log.h>

extern ChiLog chi_log;


//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void SpatialDiscretization_PWL::AddViewOfLocalContinuum(
  chi_mesh::MeshContinuum* vol_continuum,
  int num_cells,
  int* cell_indices)
{
  //================================================== Create empty view
  //                                                 for each cell
  if (!mapping_initialized)
  {
    this->cell_fe_views_mapping.reserve(vol_continuum->cells.size());
    std::vector<chi_mesh::Cell*>::iterator cellit;
    for (cellit = vol_continuum->cells.begin();
         cellit != vol_continuum->cells.end();
         cellit++)
    {
      this->cell_fe_views_mapping.push_back(-1);
    }
    mapping_initialized = true;
  }


  //================================================== Swap views for
  //                                                   specified item_id
  int cell_index = -1;
  for (int c=0; c<num_cells; c++)
  {
    cell_index = cell_indices[c];
    chi_mesh::Cell* cell = vol_continuum->cells[cell_index];

    if (cell_fe_views_mapping[cell_index]<0)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell->Type() == chi_mesh::CellType::SLAB)
      {
        auto slab_cell = dynamic_cast<chi_mesh::CellSlabV2*>(cell);
        auto cell_fe_view = new SlabFEView(slab_cell, vol_continuum);

        //cell_fe_view->PreCompute();

        this->cell_fe_views.push_back(cell_fe_view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell->Type() == chi_mesh::CellType::POLYGON)
      {
        auto poly_cell = dynamic_cast<chi_mesh::CellPolygonV2*>(cell);
        auto cell_fe_view = new PolygonFEView(poly_cell, vol_continuum, this);

        cell_fe_view->PreCompute();

        this->cell_fe_views.push_back(cell_fe_view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
      {
        auto polyh_cell = dynamic_cast<chi_mesh::CellPolyhedronV2*>(cell);
        auto cell_fe_view = new PolyhedronFEView(polyh_cell, vol_continuum, this);

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
/**Maps the cell index to a position stored locally.*/
CellFEView* SpatialDiscretization_PWL::MapFeView(int cell_glob_index)
{
  CellFEView* value = cell_fe_views.at(cell_fe_views_mapping[cell_glob_index]);
  return value;
}