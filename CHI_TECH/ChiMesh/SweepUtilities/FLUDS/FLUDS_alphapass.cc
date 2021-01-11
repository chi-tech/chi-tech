#include "FLUDS.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <ChiMesh/Cell/cell.h>

#include <chi_log.h>

extern ChiLog&     chi_log;

typedef std::vector<std::pair<int,short>> LockBox;

//###################################################################
/**Populates a flux data structure.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
InitializeAlphaElements(SPDS_ptr spds)
{
  chi_mesh::MeshContinuum*         grid = spds->grid;
  chi_mesh::sweep_management::SPLS& spls = spds->spls;

  //================================================== Initialize face
  //                                                   categorization
  num_face_categories = grid->NumberOfFaceHistogramBins();
  local_psi_stride.resize(num_face_categories,0);
  local_psi_max_elements.resize(num_face_categories,0);
  local_psi_n_block_stride.resize(num_face_categories, 0);
  local_psi_Gn_block_strideG.resize(num_face_categories,0);

  //================================================== Initialize dependent
  //                                                   locations
  size_t num_of_deplocs = spds->location_successors.size();
  deplocI_face_dof_count.resize(num_of_deplocs,0);
  deplocI_cell_views.resize(num_of_deplocs);


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
  so_cell_inco_face_face_category.reserve(spls.item_id.size());
  so_cell_outb_face_slot_indices.reserve(spls.item_id.size());
  so_cell_outb_face_face_category.reserve(spls.item_id.size());
  for (int csoi=0; csoi<spls.item_id.size(); csoi++)
  {
    int cell_local_id = spls.item_id[csoi];
    auto cell = &grid->local_cells[cell_local_id];

    local_so_cell_mapping[cell->local_id] = csoi; //Set mapping

    SlotDynamics(cell,
                 spds,
                 lock_boxes,
                 delayed_lock_box,
                 location_boundary_dependency_set);

  }//for csoi

  chi_log.Log(LOG_0VERBOSE_2) << "Done with Slot Dynamics.";
  MPI_Barrier(MPI_COMM_WORLD);



  //================================================== Populate boundary
  //                                                   dependencies
  for (auto bndry : location_boundary_dependency_set)
    boundary_dependencies.push_back(bndry);


  //                      PERFORM INCIDENT MAPPING
  //================================================== Loop over cells in
  //                                                   sweep order
  so_cell_inco_face_dof_indices.reserve(spls.item_id.size());
  for (int csoi=0; csoi<spls.item_id.size(); csoi++)
  {
    int cell_local_id = spls.item_id[csoi];
    auto cell = &grid->local_cells[cell_local_id];

    LocalIncidentMapping(cell, spds, local_so_cell_mapping);

  }//for csoi

  for (size_t fc=0; fc<num_face_categories; ++fc)
  {
    local_psi_stride[fc] = grid->GetFaceHistogramBinDOFSize(fc);
    local_psi_max_elements[fc]     = lock_boxes[fc].size();
    local_psi_n_block_stride[fc]  = local_psi_stride[fc] * lock_boxes[fc].size();
    local_psi_Gn_block_strideG[fc] = local_psi_n_block_stride[fc] * G;
  }
  delayed_local_psi_stride           = largest_face;
  delayed_local_psi_max_elements     = delayed_lock_box.size();
  delayed_local_psi_Gn_block_stride  = largest_face*delayed_lock_box.size();
  delayed_local_psi_Gn_block_strideG = delayed_local_psi_Gn_block_stride*G;

  chi_log.Log(LOG_0VERBOSE_2) << "Done with Local Incidence mapping.";
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Clean up
  so_cell_outb_face_slot_indices.shrink_to_fit();

  local_so_cell_mapping.clear();
  local_so_cell_mapping.shrink_to_fit();

  so_cell_inco_face_dof_indices.shrink_to_fit();

  nonlocal_outb_face_deplocI_slot.shrink_to_fit();

}


