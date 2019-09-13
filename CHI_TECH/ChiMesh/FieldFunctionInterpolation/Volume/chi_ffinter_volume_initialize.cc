#include "chi_ffinter_volume.h"
#include "../../Cell/cell_slab.h"
#include "../../Cell/cell_polygon.h"
#include "../../Cell/cell_polyhedron.h"


#include <chi_log.h>

extern ChiLog chi_log;

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
      if (cell->Type() == chi_mesh::CellType::SLAB)
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
      if (cell->Type() == chi_mesh::CellType::POLYGON)
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
      if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
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