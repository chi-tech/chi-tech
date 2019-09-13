#include "chi_FLUDS.h"

#include "chi_SPDS.h"

#include "../../ChiMesh/Cell/cell.h"
#include "../../ChiMesh/Cell/cell_slab.h"
#include "../../ChiMesh/Cell/cell_polygon.h"
#include "../../ChiMesh/Cell/cell_polyhedron.h"

#include <chi_log.h>

extern ChiLog     chi_log;

//###################################################################
/**Populates a flux data structure.*/
void chi_mesh::SweepManagement::FLUDS::
InitializeAlphaElements(chi_mesh::SweepManagement::SPDS* spds)
{
  chi_mesh::MeshContinuum*         grid = spds->grid;
  chi_mesh::SweepManagement::SPLS* spls = spds->spls;

  //================================================== Initialize dependent
  //                                                   locations
  int num_of_deplocs = spds->location_successors.size();
  deplocI_face_dof_count.resize(num_of_deplocs,0);
  deplocI_cell_views.resize(num_of_deplocs,std::vector<CompactCellView>());


  //                      PERFORM SLOT DYNAMICS
  //================================================== Loop over cells in
  //                                                   sweep order
  std::vector<int>  local_so_cell_mapping;
  local_so_cell_mapping.resize(grid->local_cell_glob_indices.size(),0);
  largest_face = 0;
  std::vector<std::pair<int,short>> lock_box; //cell,face index pair
  std::set<int> location_boundary_dependency_set;
  for (int csoi=0; csoi<spls->item_id.size(); csoi++)
  {
    int  cell_g_index = spls->item_id[csoi];
    auto cell         = grid->cells[cell_g_index];
    int  cell_l_index = cell->cell_local_id;

    //================================================ Create mapping of
    //                                                 local cell to sweep order
    local_so_cell_mapping[cell_l_index] = csoi;


    //================================================ Declare face_dof_mapping
    //                                                 for cell
    std::vector<std::pair<int,std::vector<short>>> inco_face_dof_mapping;

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    if (cell->Type() == chi_mesh::SLAB_CELL)
    {
      TSlab* slab_cell = (TSlab*)cell;
      SlotDynamics(slab_cell,spds,lock_box,location_boundary_dependency_set);
    }//if slab
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYGON
    else if (cell->Type() == chi_mesh::POLYGON_CELL)
    {
      TPolygon* poly_cell = (TPolygon*)cell;
      SlotDynamics(poly_cell,spds,lock_box,location_boundary_dependency_set);
    }//if polygon
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
    else if (cell->Type() == chi_mesh::POLYHEDRON_CELL)
    {
      TPolyhedron* polyh_cell = (TPolyhedron*)cell;
      SlotDynamics(polyh_cell,spds,lock_box,location_boundary_dependency_set);
    }//if polyhedron
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "InitializeAlphaElements: Unsupported cell type during "
           "first pass of function.";
      exit(EXIT_FAILURE);
    }

  }//for csoi

  //================================================== Populate boundary
  //                                                   dependencies
  std::set<int>::iterator bndry;
  for (bndry  = location_boundary_dependency_set.begin();
       bndry != location_boundary_dependency_set.end();
       bndry++)
  {
    boundary_dependencies.push_back(*bndry);
  }





  this->local_psi_stride       = largest_face;

  //                      PERFORM INCIDENT MAPPING
  //================================================== Loop over cells in
  //                                                   sweep order
  for (int csoi=0; csoi<spls->item_id.size(); csoi++)
  {
    int  cell_g_index = spls->item_id[csoi];
    auto cell         = grid->cells[cell_g_index];

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    if (cell->Type() == chi_mesh::SLAB_CELL)
    {
      TSlab* slab_cell = (TSlab*)cell;
      IncidentMapping(slab_cell,spds,local_so_cell_mapping);
    }//if slab
    // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
    else if (cell->Type() == chi_mesh::POLYGON_CELL)
    {
      TPolygon* poly_cell = (TPolygon*)cell;
      IncidentMapping(poly_cell,spds,local_so_cell_mapping);
    }//if polyhedron
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
    else if (cell->Type() == chi_mesh::POLYHEDRON_CELL)
    {
      TPolyhedron* polyh_cell = (TPolyhedron*)cell;
      IncidentMapping(polyh_cell,spds,local_so_cell_mapping);
    }//if polyhedron

  }//for csoi

  this->local_psi_stride       = largest_face;
  this->local_psi_max_elements = lock_box.size();
  this->local_psi_Gn_block_stride = largest_face*lock_box.size();
  this->local_psi_Gn_block_strideG = local_psi_Gn_block_stride*G;



  //================================================== Clean up
  this->so_cell_outb_face_slot_indices.shrink_to_fit();

  local_so_cell_mapping.clear();
  local_so_cell_mapping.shrink_to_fit();

  so_cell_inco_face_dof_indices.shrink_to_fit();

  nonlocal_outb_face_deplocI_slot.shrink_to_fit();

}