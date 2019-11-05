#include "pwl.h"

#include "../../../ChiMesh/Cell/cell_triangle.h"
#include "CellViews/pwl_slab.h"
#include "CellViews/pwl_triangle.h"
#include "CellViews/pwl_polygon.h"
#include "CellViews/pwl_polyhedron.h"
#include<typeinfo>

#include <chi_log.h>
#include <chi_mpi.h>

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
      //========================================= If slab item_id
      if (cell->Type() == chi_mesh::CellType::SLAB)
      {
        SlabFEView* view =
          new SlabFEView((chi_mesh::CellSlab*)cell,vol_continuum);

        this->cell_fe_views.push_back(view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }

      //========================================= If triangle item_id
      else if (typeid(*(cell)) == typeid(chi_mesh::CellTriangle) )
      {
        TriangleFEView* view =
          new TriangleFEView((chi_mesh::CellTriangle*)(cell),vol_continuum);

        this->cell_fe_views.push_back(view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }

      //========================================= If polygon item_id
      else if (cell->Type() == chi_mesh::CellType::POLYGON)
      {
        PolygonFEView* view =
          new PolygonFEView((chi_mesh::CellPolygon*)(cell),vol_continuum,this);

        view->PreCompute();
        this->cell_fe_views.push_back(view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }

      //========================================= If polyhedron item_id
      else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
      {
        PolyhedronFEView* view =
          new PolyhedronFEView(
            (chi_mesh::CellPolyhedron*)(cell),
            vol_continuum,
            this);

        view->PreCompute();
        view->CleanUp();
        this->cell_fe_views.push_back(view);
        cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
      }
        //========================================= If new_base
      else if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
      {
        auto cell_base = dynamic_cast<chi_mesh::CellBase*>(cell);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLABV2
        if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
        {
          auto slab_cell = dynamic_cast<chi_mesh::CellSlabV2*>(cell_base);
          auto cell_fe_view = new SlabFEView(slab_cell, vol_continuum);

          //cell_fe_view->PreCompute();

          this->cell_fe_views.push_back(cell_fe_view);
          cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGONV2
        else if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
        {
          auto poly_cell = dynamic_cast<chi_mesh::CellPolygonV2*>(cell_base);
          auto cell_fe_view = new PolygonFEView(poly_cell, vol_continuum, this);

          cell_fe_view->PreCompute();

          this->cell_fe_views.push_back(cell_fe_view);
          cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRONV2
        else if (cell_base->Type2() == chi_mesh::CellType::POLYHEDRONV2)
        {
          auto polyh_cell = dynamic_cast<chi_mesh::CellPolyhedronV2*>(cell_base);
          auto cell_fe_view = new PolyhedronFEView(polyh_cell, vol_continuum, this);

          cell_fe_view->PreCompute();
          cell_fe_view->CleanUp();
          this->cell_fe_views.push_back(cell_fe_view);
          cell_fe_views_mapping[cell_index] = this->cell_fe_views.size()-1;
        }
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