#ifndef AUX_FLUDS_h
#define AUX_FLUDS_h

#include "FLUDS.h"

//########################################################################
/**Auxiliary flux data structure.*/
class chi_mesh::sweep_management::AUX_FLUDS :
  public chi_mesh::sweep_management::FLUDS
{
public:
  AUX_FLUDS(PRIMARY_FLUDS& set_primary, int in_G);

  //Inherited from base FLUDS
//public:
//  std::vector<size_t>& local_psi_stride;               //face category face size
//  std::vector<size_t>& local_psi_max_elements;         //number of face in each cat
//  const size_t         delayed_local_psi_stride;       //face category face size
//  const size_t         delayed_local_psi_max_elements; //number of face in each cat
//  const size_t         num_face_categories;


private:
  const int largest_face;
  const int G;

  //local_psi_Gn_block_stride[fc]. Given face category fc, the value is
  //total number of faces that store information in this categorie's buffer
  std::vector<size_t>& local_psi_Gn_block_stride;
  std::vector<size_t>  local_psi_Gn_block_strideG;
  const size_t delayed_local_psi_Gn_block_stride;
  size_t delayed_local_psi_Gn_block_strideG;

  //======================================== References to psi vectors
  //ref_local_psi[fc]. Each category fc has its own
  //  local psi interface vector storing all aggregated angles and groups.
  //ref_delayed_local_psi. All delayed information is in this vector.
  //ref_deplocI_outgoing_psi[deplocI]. Each dependent location I has
  //  its own interface vector
  //ref_prelocI_outgoing_psi[prelocI]. Each predecessor location I has
  //  its own interface vector
  //ref_delayed_prelocI_outgoing_psi[prelocI]. Each delayed predecessor
  //  location I has its own interface vector
  std::vector<std::vector<double>>*  ref_local_psi;
  std::vector<double>*               ref_delayed_local_psi;
  std::vector<double>*               ref_delayed_local_psi_old;
  std::vector<std::vector<double>>*  ref_deplocI_outgoing_psi;
  std::vector<std::vector<double>>*  ref_prelocI_outgoing_psi;
  std::vector<std::vector<double>>*  ref_boundryI_incoming_psi;

  std::vector<std::vector<double>>*  ref_delayed_prelocI_outgoing_psi;
  std::vector<std::vector<double>>*  ref_delayed_prelocI_outgoing_psi_old;
private:
  //======================================== Alpha elements

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the slot address in the local psi vector where the first
  // face dof will store its data
  std::vector<std::vector<int>>& so_cell_outb_face_slot_indices;

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<std::vector<int>>& so_cell_outb_face_face_category;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // that will hold a pair. Pair-first holds the slot address where this
  // face's upwind data is stored. Pair-second is a mapping of
  // each of this face's dofs to the upwinded face's dofs
  std::vector<std::vector<std::pair<int,std::vector<int>> >>&
    so_cell_inco_face_dof_indices;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<std::vector<int>>& so_cell_inco_face_face_category;

  //Inherited from base FLUDS
//public:
//  // This is a small vector [deplocI] that holds the number of
//  // face dofs for each dependent location.
//  std::vector<int>&    deplocI_face_dof_count;
//
//  // Very small vector listing the boundaries this location depends on
//  std::vector<int>& boundary_dependencies;

private:
  // This is a vector [non_local_outgoing_face_count]
  // that maps a face to a dependent location and associated slot index
  std::vector<std::pair<int,int>>& nonlocal_outb_face_deplocI_slot;

  //======================================== Beta elements
  //Inherited from base FLUDS
//public:
//  // This is a small vector [prelocI] that holds the number of
//  // face dofs for each predecessor location.
//  std::vector<int>&    prelocI_face_dof_count;
//
//  std::vector<int>&    delayed_prelocI_face_dof_count;


private:
  // This is a vector [nonlocal_inc_face_counter] containing
  // AlphaPairs. AlphaPair-first is the prelocI index and
  // AlphaPair-second is a BetaPair. The BetaPair-first is the slot where
  // the face storage begins and BetaPair-second is a dof mapping
  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>&
    nonlocal_inc_face_prelocI_slot_dof;

  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>&
    delayed_nonlocal_inc_face_prelocI_slot_dof;

public:
  /**Passes pointers from sweep buffers to FLUDS so
   * that chunk utilities function as required. */
  void SetReferencePsi(
    std::vector<std::vector<double>>*  local_psi,
    std::vector<double>*               delayed_local_psi,
    std::vector<double>*               delayed_local_psi_old,
    std::vector<std::vector<double>>*  deplocI_outgoing_psi,
    std::vector<std::vector<double>>*  prelocI_outgoing_psi,
    std::vector<std::vector<double>>*  boundryI_incoming_psi,
    std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi,
    std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi_old)
  override
  {
    ref_local_psi = local_psi;
    ref_delayed_local_psi = delayed_local_psi;
    ref_delayed_local_psi_old = delayed_local_psi_old;
    ref_deplocI_outgoing_psi = deplocI_outgoing_psi;
    ref_prelocI_outgoing_psi = prelocI_outgoing_psi;
    ref_boundryI_incoming_psi = boundryI_incoming_psi;
    ref_delayed_prelocI_outgoing_psi = delayed_prelocI_outgoing_psi;
    ref_delayed_prelocI_outgoing_psi_old = delayed_prelocI_outgoing_psi_old;
  }

  double*  OutgoingPsi(int cell_so_index, int outb_face_counter,
                       int face_dof, int n) override;
  double*  UpwindPsi(int cell_so_index, int inc_face_counter,
                     int face_dof,int g, int n) override;


  double*  NLOutgoingPsi(int outb_face_count,int face_dof, int n) override;

  double*  NLUpwindPsi(int nonl_inc_face_counter,
                       int face_dof,int g, int n) override;
};

#endif