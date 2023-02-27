#ifndef CHI_FLUDS_H
#define CHI_FLUDS_H

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

namespace chi_mesh
{
  class GridFaceHistogram;
}

#include <iostream>

//face_slot index, vertex ids
typedef std::pair<int,std::vector<uint64_t>>             CompactFaceView;

//cell_global_id, faces
typedef std::pair<int,std::vector<CompactFaceView>> CompactCellView;

namespace chi_mesh::sweep_management
{
  struct FaceNodalMapping
  {
    const int associated_face;
    const std::vector<short> node_mapping;

    FaceNodalMapping(int in_ass_face,
                     std::vector<short>& in_node_mapping) :
                     associated_face(in_ass_face),
                     node_mapping(in_node_mapping) {}
  };
  typedef std::vector<FaceNodalMapping> CellFaceNodalMapping;

  class FLUDS
  {
  public:
    size_t              num_face_categories=0;            //Number of face categories
    std::vector<size_t> local_psi_stride;                 //Group-angle-faceDOF stride per cat
    std::vector<size_t> local_psi_max_elements;           //Number of faces in each cat
    size_t              delayed_local_psi_stride=0;       //Group-angle-faceDOF stride delayed cat
    size_t              delayed_local_psi_max_elements=0; //Number of faces in delayed cat

    std::shared_ptr<SPDS> spds_;

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
    virtual
    void SetReferencePsi(
      std::vector<std::vector<double>>*  local_psi,
      std::vector<double>*               delayed_local_psi,
      std::vector<double>*               delayed_local_psi_old,
      std::vector<std::vector<double>>*  deplocI_outgoing_psi,
      std::vector<std::vector<double>>*  prelocI_outgoing_psi,
      std::vector<std::vector<double>>*  boundryI_incoming_psi,
      std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi,
      std::vector<std::vector<double>>*  delayed_prelocI_outgoing_psi_old)=0;

    virtual
    double*  OutgoingPsi(int cell_so_index, int outb_face_counter,
                         int face_dof, int n) = 0;
    virtual
    double*  UpwindPsi(int cell_so_index, int inc_face_counter,
                       int face_dof,int g, int n) = 0;

    virtual
    double*  NLOutgoingPsi(int nonl_outb_face_counter, int face_dof, int n) = 0;

    virtual
    double*  NLUpwindPsi(int nonl_inc_face_counter,
                         int face_dof,int g, int n) = 0;

    virtual ~FLUDS()=default;
  };

  struct INCOMING_FACE_INFO
  {
    int slot_address=0;
    short* upwind_dof_mapping = nullptr;

    void Setup(const std::pair<int,std::vector<short>>& input)
    {
      slot_address = input.first;
      upwind_dof_mapping = new short[input.second.size()];
      std::copy(input.second.begin(),input.second.end(),upwind_dof_mapping);
    }

    ~INCOMING_FACE_INFO()
    {
      delete [] upwind_dof_mapping;
    }
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
//  size_t              num_face_categories=0;            //Number of face categories
//  std::vector<size_t> local_psi_stride;                 //Group-angle-faceDOF stride per cat
//  std::vector<size_t> local_psi_max_elements;           //Number of faces in each cat
//  size_t              delayed_local_psi_stride=0;       //Group-angle-faceDOF stride delayed cat
//  size_t              delayed_local_psi_max_elements=0; //Number of faces in delayed cat
//
//public:
//  // This is a small vector [deplocI] that holds the number of
//  // face dofs for each dependent location.
//  std::vector<int>    deplocI_face_dof_count;
//
//  // Very small vector listing the boundaries this location depends on
//  std::vector<int> boundary_dependencies;
//
//public:
//  // This is a small vector [prelocI] that holds the number of
//  // face dofs for each predecessor location.
//  std::vector<int>    prelocI_face_dof_count;
//
//  std::vector<int>    delayed_prelocI_face_dof_count;
protected:
  const SPDS& spds_;
  const GridFaceHistogram& grid_face_histogram_;

private:
  int largest_face=0;
  size_t G=0;

  //local_psi_n_block_stride[fc]. Given face category fc, the value is
  //total number of faces that store information in this category's buffer
  //per angle
  std::vector<size_t> local_psi_n_block_stride;
  std::vector<size_t> local_psi_Gn_block_strideG;
  size_t delayed_local_psi_Gn_block_stride=0;
  size_t delayed_local_psi_Gn_block_strideG=0;

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
  std::vector<std::vector<double>>*  ref_local_psi = nullptr;
  std::vector<double>*               ref_delayed_local_psi = nullptr;
  std::vector<double>*               ref_delayed_local_psi_old = nullptr;
  std::vector<std::vector<double>>*  ref_deplocI_outgoing_psi = nullptr;
  std::vector<std::vector<double>>*  ref_prelocI_outgoing_psi = nullptr;
  std::vector<std::vector<double>>*  ref_boundryI_incoming_psi = nullptr;

  std::vector<std::vector<double>>*  ref_delayed_prelocI_outgoing_psi = nullptr;
  std::vector<std::vector<double>>*  ref_delayed_prelocI_outgoing_psi_old = nullptr;
private:
  //======================================== Alpha elements

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the slot address in the local psi vector where the first
  // face dof will store its data
  std::vector<int*>
    so_cell_outb_face_slot_indices;

  // This is a vector [cell_sweep_order_index][outgoing_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<short*>
    so_cell_outb_face_face_category;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // that will hold a structure. struct.slot_address holds the slot address
  // where this face's upwind data is stored. struct.upwind_dof_mapping is
  // a mapping of each of this face's dofs to the upwinded face's dofs
private:
  std::vector<INCOMING_FACE_INFO*>
    so_cell_inco_face_dof_indices;

  // This is a vector [cell_sweep_order_index][incoming_face_count]
  // which holds the face categorization for the face. i.e. the local
  // psi vector that hold faces of the same category.
  std::vector<short*>
    so_cell_inco_face_face_category;

private:
  // This is a vector [non_local_outgoing_face_count]
  // that maps a face to a dependent location and associated slot index
  std::vector<std::pair<int,int>>
    nonlocal_outb_face_deplocI_slot;

  //======================================== Beta elements
private:
  // This is a vector [nonlocal_inc_face_counter] containing
  // AlphaPairs. AlphaPair-first is the prelocI index and
  // AlphaPair-second is a BetaPair. The BetaPair-first is the slot where
  // the face storage begins and BetaPair-second is a dof mapping
  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>
    nonlocal_inc_face_prelocI_slot_dof;

  std::vector<std::pair<int,std::pair<int,std::vector<int>>>>
    delayed_nonlocal_inc_face_prelocI_slot_dof;

  //======================================== Temporary Utility Members
private:
  // This is a vector [dependent_location][unordered_cell_index]
  // that holds an AlphaPair. AlphaPair-first is the cell's global_id
  // and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  // comprises BetaPair-first = face_slot_index (the location of this
  // faces data in the psi vector, and then a vector of vertex indexes
  // that can be used for dof_mapping.
  // Filled during slot-dynamics.
  // Cleared after beta-pass.
  std::vector<std::vector<CompactCellView>> deplocI_cell_views;

private:
  // This is a vector [predecessor_location][unordered_cell_index]
  // that holds an AlphaPair. AlphaPair-first is the cell's global_id
  // and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  // comprises BetaPair-first = face_slot_index (the location of this
  // faces data in the psi vector, and then a vector of vertex indexes
  // that can be used for dof_mapping.
  // Filled in beta-pass
  // Cleared after beta-pass.
  std::vector<std::vector<CompactCellView>> prelocI_cell_views;
  std::vector<std::vector<CompactCellView>> delayed_prelocI_cell_views;

  std::vector<CellFaceNodalMapping>& grid_nodal_mappings;

public:
  PRIMARY_FLUDS(size_t in_G,
                std::vector<CellFaceNodalMapping>& in_grid_nodal_mappings,
                const SPDS& spds,
                const chi_mesh::GridFaceHistogram& grid_face_histogram);

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

public:
  typedef std::shared_ptr<chi_mesh::sweep_management::SPDS> SPDS_ptr;
  //alphapass.cc
  void InitializeAlphaElements(const SPDS& spds);

  void AddFaceViewToDepLocI(int deplocI, int cell_g_index,
                            int face_slot, const chi_mesh::CellFace& face);

  //alphapass_slotdynamics.cc
  void SlotDynamics(const chi_mesh::Cell& cell,
                    const SPDS& spds,
                    std::vector<std::vector<std::pair<int,short>>>& lock_boxes,
                    std::vector<std::pair<int,short>>& delayed_lock_box,
                    std::set<int>& location_boundary_dependency_set);
  //alphapass_inc_mapping.cc
  void LocalIncidentMapping(const chi_mesh::Cell& cell,
                            const SPDS& spds,
                            std::vector<int>&  local_so_cell_mapping);

  //betapass.cc
  void InitializeBetaElements(const SPDS& spds,
                              int tag_index=0);

  void SerializeCellInfo(std::vector<CompactCellView>* cell_views,
                         std::vector<int>& face_indices,
                         int num_face_dofs);
  void DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                           std::vector<int>* face_indices,
                           int& num_face_dofs);

  //betapass_nonlocal_inc_mapping.cc
  void NonLocalIncidentMapping(const chi_mesh::Cell& cell,
                               const SPDS& spds);

  //FLUDS_chunk_utilities.cc
  double*  OutgoingPsi(int cell_so_index, int outb_face_counter,
                       int face_dof, int n) override;
  double*  UpwindPsi(int cell_so_index, int inc_face_counter,
                     int face_dof,int g, int n) override;


  double*  NLOutgoingPsi(int outb_face_count,int face_dof, int n) override;

  double*  NLUpwindPsi(int nonl_inc_face_counter,
                       int face_dof,int g, int n) override;

  ~PRIMARY_FLUDS() override
  {
    for (auto& val : so_cell_outb_face_slot_indices) delete [] val;
    for (auto& val : so_cell_outb_face_face_category) delete [] val;
    for (auto& val : so_cell_inco_face_dof_indices) delete [] val;
    for (auto& val : so_cell_inco_face_face_category) delete [] val;

  }

};

#endif //CHI_FLUDS_H