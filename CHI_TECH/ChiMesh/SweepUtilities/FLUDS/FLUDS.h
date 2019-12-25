#ifndef _chi_FLUDS_h
#define _chi_FLUDS_h


#include <stdio.h>

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include <ChiMesh/Cell/cell.h>

//face_slot index, vertex ids
typedef std::pair<int,std::vector<int>>             CompactFaceView;

//cell_global_id, faces
typedef std::pair<int,std::vector<CompactFaceView>> CompactCellView;

namespace chi_mesh::sweep_management
{
  class FLUDS
  {
  public:
    std::vector<size_t> local_psi_stride;               //face category face size
    std::vector<size_t> local_psi_max_elements;         //number of face in each cat
    size_t              delayed_local_psi_stride;       //face category face size
    size_t              delayed_local_psi_max_elements; //number of face in each cat
    size_t              num_face_categories;

  public:
    // This is a small vector [deplocI] that holds the number of
    // face dofs for each dependent location.
    std::vector<int>    deplocI_face_dof_count;

    // Very small vector listing the boundaries this location depends on
    std::vector<int> boundary_dependencies;

  public:
    // This is a small vector [prelocI] that holds the number of
    // face dofs for each predecessor location.
    std::vector<int>    prelocI_face_dof_count;

    std::vector<int>    delayed_prelocI_face_dof_count;


  public:
    virtual ~FLUDS() {};
    virtual
    void SetReferencePsi(
      std::vector<std::vector<double>>*  local_psi,
       std::vector<double>*               delayed_local_psi,
       std::vector<std::vector<double>>*  deplocI_outgoing_psi,
       std::vector<std::vector<double>>*  prelocI_outgoing_psi,
       std::vector<std::vector<double>>*  boundryI_incoming_psi,
       std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi)=0;

    virtual
    double*  OutgoingPsi(int cell_so_index, int outb_face_counter,
                         int face_dof, int n) = 0;
    virtual
    double*  UpwindPsi(int cell_so_index, int inc_face_counter,
                       int face_dof,int g, int n) = 0;

    virtual
    double*  NLOutgoingPsi(int outb_face_count,int face_dof, int n) = 0;

    virtual
    double*  NLUpwindPsi(int nonl_inc_face_counter,
                         int face_dof,int g, int n) = 0;
  };
}

//###################################################################
/**Improved Flux Data Structure (FLUDS).*/
class chi_mesh::sweep_management::PRIMARY_FLUDS :
  public chi_mesh::sweep_management::FLUDS
{
  /**\relates chi_mesh::sweep_management::AUX_FLUDS*/
  friend chi_mesh::sweep_management::AUX_FLUDS;

  //Inherited from base FLUDS
//public:
//  std::vector<size_t> local_psi_stride;               //face category face size
//  std::vector<size_t> local_psi_max_elements;         //number of face in each cat
//  size_t              delayed_local_psi_stride;       //face category face size
//  size_t              delayed_local_psi_max_elements; //number of face in each cat
//  size_t              num_face_categories;


private:
  int largest_face;
  int G;

  //local_psi_Gn_block_stride[fc]. Given face category fc, the value is
  //total number of faces that store information in this categorie's buffer
  std::vector<size_t> local_psi_Gn_block_stride;
  std::vector<size_t> local_psi_Gn_block_strideG;
  size_t delayed_local_psi_Gn_block_stride;
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
  std::vector<std::vector<double>>*  ref_deplocI_outgoing_psi;
  std::vector<std::vector<double>>*  ref_prelocI_outgoing_psi;
  std::vector<std::vector<double>>*  ref_boundryI_incoming_psi;

  std::vector<std::vector<double>>*  ref_delayed_prelocI_outgoing_psi;
private:
  //======================================== Alpha elements

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the slot address in the local psi vector where the first
  // face dof will store its data
  std::vector<std::vector<int>>
    so_cell_outb_face_slot_indices;

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<std::vector<int>>
    so_cell_outb_face_face_category;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // that will hold a pair. Pair-first holds the slot address where this
  // face's upwind data is stored. Pair-second is a mapping of
  // each of this face's dofs to the upwinded face's dofs
  std::vector<std::vector<std::pair<int,std::vector<int>> >>
    so_cell_inco_face_dof_indices;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<std::vector<int>>
    so_cell_inco_face_face_category;

  //Inherited from base FLUDS
//public:
//  // This is a small vector [deplocI] that holds the number of
//  // face dofs for each dependent location.
//  std::vector<int>    deplocI_face_dof_count;
//
//  // Very small vector listing the boundaries this location depends on
//  std::vector<int> boundary_dependencies;

private:
  // This is a vector [non_local_outgoing_face_count]
  // that maps a face to a dependent location and associated slot index
  std::vector<std::pair<int,int>>
    nonlocal_outb_face_deplocI_slot;

private:
  // This is a vector [dependent_location][unordered_cell_index]
  // that holds an AlphaPair. AlphaPair-first is the cell's global_id
  // and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  // comprises BetaPair-first = face_slot_index (the location of this
  // faces data in the psi vector, and then a vector of vertex indexes
  // that can be used for dof_mapping.
  std::vector<std::vector<CompactCellView>>
    deplocI_cell_views;

  //======================================== Beta elements
  //Inherited from base FLUDS
//public:
//  // This is a small vector [prelocI] that holds the number of
//  // face dofs for each predecessor location.
//  std::vector<int>    prelocI_face_dof_count;
//
//  std::vector<int>    delayed_prelocI_face_dof_count;

private:
  // This is a vector [predecessor_location][unordered_cell_index]
  // that holds an AlphaPair. AlphaPair-first is the cell's global_id
  // and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  // comprises BetaPair-first = face_slot_index (the location of this
  // faces data in the psi vector, and then a vector of vertex indexes
  // that can be used for dof_mapping.
  std::vector<std::vector<CompactCellView>>
    prelocI_cell_views;

  std::vector<std::vector<CompactCellView>>
      delayed_prelocI_cell_views;

private:
  // This is a vector [nonlocal_inc_face_counter] containing
  // AlphaPairs. AlphaPair-first is the prelocI index and
  // AlphaPair-second is a BetaPair. The BetaPair-first is the slot where
  // the face storage begins and BetaPair-second is a dof mapping
  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>
    nonlocal_inc_face_prelocI_slot_dof;

  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>
      delayed_nonlocal_inc_face_prelocI_slot_dof;

public:
  PRIMARY_FLUDS(int in_G): G(in_G)
  { }

public:
  /**Passes pointers from sweep buffers to FLUDS so
   * that chunk utilities function as required. */
  void SetReferencePsi(
    std::vector<std::vector<double>>*  local_psi,
    std::vector<double>*               delayed_local_psi,
    std::vector<std::vector<double>>*  deplocI_outgoing_psi,
    std::vector<std::vector<double>>*  prelocI_outgoing_psi,
    std::vector<std::vector<double>>*  boundryI_incoming_psi,
    std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi)
    override
  {
    ref_local_psi = local_psi;
    ref_delayed_local_psi = delayed_local_psi;
    ref_deplocI_outgoing_psi = deplocI_outgoing_psi;
    ref_prelocI_outgoing_psi = prelocI_outgoing_psi;
    ref_boundryI_incoming_psi = boundryI_incoming_psi;
    ref_delayed_prelocI_outgoing_psi = delayed_prelocI_outgoing_psi;
  }

public:
  //alphapass.cc
  void InitializeAlphaElements(chi_mesh::sweep_management::SPDS *spds);

  void AddFaceViewToDepLocI(int deplocI, int cell_g_index,
                            int face_slot, chi_mesh::CellFace& face);

  //alphapass_slotdynamics.cc
  void SlotDynamics(chi_mesh::Cell *cell,
                    chi_mesh::sweep_management::SPDS* spds,
                    std::vector<std::vector<std::pair<int,short>>>& lock_boxes,
                    std::vector<std::pair<int,short>>& delayed_lock_box,
                    std::set<int>& location_boundary_dependency_set);
  //alphapass_inc_mapping.cc
  void LocalIncidentMapping(chi_mesh::Cell *cell,
                            chi_mesh::sweep_management::SPDS* spds,
                            std::vector<int>&  local_so_cell_mapping);

  //betapass.cc
  void InitializeBetaElements(chi_mesh::sweep_management::SPDS *spds,
                              int tag_index=0);

  void SerializeCellInfo(std::vector<CompactCellView>* cell_views,
                         std::vector<int>& face_indices,
                         int num_face_dofs);
  void DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                           std::vector<int>* face_indices,
                           int& num_face_dofs);

  //betapass_nonlocal_inc_mapping.cc
  void NonLocalIncidentMapping(chi_mesh::Cell *cell,
                               chi_mesh::sweep_management::SPDS* spds);

  //FLUDS_chunk_utilities.cc
  double*  OutgoingPsi(int cell_so_index, int outb_face_counter,
                       int face_dof, int n) override;
  double*  UpwindPsi(int cell_so_index, int inc_face_counter,
                     int face_dof,int g, int n) override;


  double*  NLOutgoingPsi(int outb_face_count,int face_dof, int n) override;

  double*  NLUpwindPsi(int nonl_inc_face_counter,
                       int face_dof,int g, int n) override;

};

#endif