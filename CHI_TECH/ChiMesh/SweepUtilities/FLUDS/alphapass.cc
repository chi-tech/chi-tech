#include "FLUDS.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <ChiMesh/Cell/cell.h>

#include <chi_log.h>

extern ChiLog     chi_log;

typedef std::vector<std::pair<int,short>> LockBox;

//###################################################################
/**Populates a flux data structure.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
InitializeAlphaElements(chi_mesh::sweep_management::SPDS* spds)
{
  chi_mesh::MeshContinuum*         grid = spds->grid;
  chi_mesh::sweep_management::SPLS* spls = spds->spls;

  num_face_categories = grid->NumberOfFaceHistogramBins();
  local_psi_stride.resize(num_face_categories,0);
  local_psi_max_elements.resize(num_face_categories,0);
  local_psi_Gn_block_stride.resize(num_face_categories,0);
  local_psi_Gn_block_strideG.resize(num_face_categories,0);


  //================================================== Initialize dependent
  //                                                   locations
  size_t num_of_deplocs = spds->location_successors.size();
  deplocI_face_dof_count.resize(num_of_deplocs,0);
  deplocI_cell_views.resize(num_of_deplocs,std::vector<CompactCellView>());


  //                      PERFORM SLOT DYNAMICS
  //================================================== Loop over cells in
  //                                                   sweep order
  // Given a local cell index, gives the so index
  std::vector<int>  local_so_cell_mapping;
  local_so_cell_mapping.resize(grid->local_cell_glob_indices.size(),0);

  largest_face = 0; // Will contain the max dofs per face
  std::vector<LockBox> lock_boxes(num_face_categories); //cell,face index pairs
  LockBox              delayed_lock_box;
  std::set<int> location_boundary_dependency_set;

  // csoi = cell sweep order index
  for (int csoi=0; csoi<spls->item_id.size(); csoi++)
  {
    int  cell_g_index = spls->item_id[csoi];           // Global index
    auto cell         = grid->cells[cell_g_index];
    local_so_cell_mapping[cell->cell_local_id] = csoi; //Set mapping

    SlotDynamics(cell,spds,lock_boxes,delayed_lock_box,location_boundary_dependency_set);

  }//for csoi


  //================================================== Populate boundary
  //                                                   dependencies
  for (auto bndry : location_boundary_dependency_set)
    boundary_dependencies.push_back(bndry);

  //                      PERFORM INCIDENT MAPPING
  //================================================== Loop over cells in
  //                                                   sweep order
  for (int csoi=0; csoi<spls->item_id.size(); csoi++)
  {
    int  cell_g_index = spls->item_id[csoi];
    auto cell         = grid->cells[cell_g_index];

    LocalIncidentMapping(cell, spds, local_so_cell_mapping);

  }//for csoi

  for (size_t fc=0; fc<num_face_categories; ++fc)
  {
    local_psi_stride[fc] = grid->GetFaceHistogramBinDOFSize(fc);
    local_psi_max_elements[fc]     = lock_boxes[fc].size();
    local_psi_Gn_block_stride[fc]  = local_psi_stride[fc]*lock_boxes[fc].size();
    local_psi_Gn_block_strideG[fc] = local_psi_Gn_block_stride[fc]*G;
  }
  delayed_local_psi_stride           = largest_face;
  delayed_local_psi_max_elements     = delayed_lock_box.size();
  delayed_local_psi_Gn_block_stride  = largest_face*delayed_lock_box.size();
  delayed_local_psi_Gn_block_strideG = delayed_local_psi_Gn_block_stride*G;


  //================================================== Clean up
  this->so_cell_outb_face_slot_indices.shrink_to_fit();

  local_so_cell_mapping.clear();
  local_so_cell_mapping.shrink_to_fit();

  so_cell_inco_face_dof_indices.shrink_to_fit();

  nonlocal_outb_face_deplocI_slot.shrink_to_fit();

}

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector.*/
void  chi_mesh::sweep_management::PRIMARY_FLUDS::
AddFaceViewToDepLocI(int deplocI, int cell_g_index, int face_slot,
                     chi_mesh::CellFace& face)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (int c=0; c<deplocI_cell_views[deplocI].size(); c++)
  {
    if (deplocI_cell_views[deplocI][c].first == cell_g_index)
    {
      cell_already_there = true;
      deplocI_cell_views[deplocI][c].second.
        emplace_back(face_slot, face.vertex_ids);
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.
      emplace_back(face_slot, face.vertex_ids);

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }


}
