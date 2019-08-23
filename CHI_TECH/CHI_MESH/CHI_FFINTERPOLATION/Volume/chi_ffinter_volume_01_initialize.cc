#include "chi_ffinter_volume.h"
#include "../../CHI_CELL/cell_slab.h"
#include "../../CHI_CELL/cell_polygon.h"
#include "../../CHI_CELL/cell_polyhedron.h"


#include <chi_log.h>

extern CHI_LOG chi_log;

//###################################################################
/**Initializes the volume field function interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Initialize()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing volume interpolator.";
  //================================================== Check grid available
  if (field_functions.size() == 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unassigned field function in volume field function interpolator.";
    exit(EXIT_FAILURE);
  } else
  {
    this->grid_view = field_functions[0]->grid;
  }

  //================================================== Find cell inside volume
  size_t num_local_cells = grid_view->local_cell_glob_indices.size();
  for (int lc=0; lc<num_local_cells; lc++)
  {
    int cell_glob_index = grid_view->local_cell_glob_indices[lc];
    auto cell = grid_view->cells[cell_glob_index];

    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell->centroid);

    if (inside_logvolume)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slab
      if (typeid(*cell) == typeid(chi_mesh::CellSlab))
      {
        chi_mesh::CellSlab* slab_cell = (chi_mesh::CellSlab*)cell;

        for (int i=0; i<2; i++)
        {
          cfem_local_nodes_needed_unmapped.push_back(slab_cell->v_indices[i]);
          pwld_local_nodes_needed_unmapped.push_back(i);
          pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
        }//for dof
      }

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
      {
        chi_mesh::CellPolygon* poly_cell = (chi_mesh::CellPolygon*)cell;

        for (int i=0; i<poly_cell->v_indices.size(); i++)
        {
          cfem_local_nodes_needed_unmapped.push_back(poly_cell->v_indices[i]);
          pwld_local_nodes_needed_unmapped.push_back(i);
          pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
        }//for dof
      }

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
      {
        chi_mesh::CellPolyhedron* polyh_cell = (chi_mesh::CellPolyhedron*)cell;

        for (int i=0; i<polyh_cell->v_indices.size(); i++)
        {
          cfem_local_nodes_needed_unmapped.push_back(polyh_cell->v_indices[i]);
          pwld_local_nodes_needed_unmapped.push_back(i);
          pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
        }//for dof
      }//if Polyhedron
    }//if inside logicalVol

  }//for local cell
}