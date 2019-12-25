#include "AUX_FLUDS.h"

#include <chi_log.h>
extern ChiLog     chi_log;


//######################################################### Constr
/**This constructor initializes an auxiliary FLUDS based
 * on a primary FLUDS. The restriction here is that the
 * auxiliary FLUDS has the exact same sweep ordering as the
 * primary FLUDS.*/
chi_mesh::sweep_management::AUX_FLUDS::
  AUX_FLUDS(PRIMARY_FLUDS &primary, int in_G) :
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Initializing references
  largest_face( primary.largest_face ),
  G( in_G ),

  local_psi_Gn_block_stride( primary.local_psi_Gn_block_stride ),

  delayed_local_psi_Gn_block_stride( primary.delayed_local_psi_Gn_block_stride ),
  //================ Alpha Elements
  so_cell_outb_face_slot_indices( primary.so_cell_outb_face_slot_indices ),
  so_cell_outb_face_face_category( primary.so_cell_outb_face_face_category ),

  so_cell_inco_face_dof_indices( primary.so_cell_inco_face_dof_indices ),
  so_cell_inco_face_face_category( primary.so_cell_inco_face_face_category ),

  //================ Beta Elements
  nonlocal_outb_face_deplocI_slot( primary.nonlocal_outb_face_deplocI_slot ),

  nonlocal_inc_face_prelocI_slot_dof( primary.nonlocal_inc_face_prelocI_slot_dof ),
  delayed_nonlocal_inc_face_prelocI_slot_dof(
    primary.delayed_nonlocal_inc_face_prelocI_slot_dof )

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Actual constructor
{
  //============================== Small copied items
  //                               defined on base class
  local_psi_stride               = primary.local_psi_stride;
  local_psi_max_elements         = primary.local_psi_max_elements;
  delayed_local_psi_stride       = primary.delayed_local_psi_stride;
  delayed_local_psi_max_elements = primary.delayed_local_psi_max_elements;
  num_face_categories            = primary.num_face_categories;

  deplocI_face_dof_count         = primary.deplocI_face_dof_count;
  boundary_dependencies          = primary.boundary_dependencies;

  prelocI_face_dof_count         = primary.prelocI_face_dof_count;
  delayed_prelocI_face_dof_count = primary.delayed_prelocI_face_dof_count;

  //============================== Adjusting for different group aggregate
  for (auto& val : local_psi_Gn_block_stride)
    local_psi_Gn_block_strideG.push_back(val*G);

  delayed_local_psi_Gn_block_strideG = delayed_local_psi_Gn_block_stride*G;

}



//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector
 * and returns a reference to it.*/
double*  chi_mesh::sweep_management::AUX_FLUDS::
OutgoingPsi(int cell_so_index, int outb_face_counter,
            int face_dof, int n)
{
  // Face category
  int fc = so_cell_outb_face_face_category[cell_so_index][outb_face_counter];

  if (fc >= 0)
  {
    size_t index =
      local_psi_Gn_block_strideG[fc]*n +
      so_cell_outb_face_slot_indices[cell_so_index][outb_face_counter]*
      local_psi_stride[fc]*G +
      face_dof*G;

    return &(ref_local_psi->operator[](fc))[index];
  }
  else
  {
    size_t index =
      delayed_local_psi_Gn_block_strideG*n +
      so_cell_outb_face_slot_indices[cell_so_index][outb_face_counter]*
      delayed_local_psi_stride*G +
      face_dof*G;

    return &ref_delayed_local_psi->operator[](index);
  }

}

//###################################################################
/**Given a outbound face counter this method returns a pointer
 * to the location*/
double*  chi_mesh::sweep_management::AUX_FLUDS::
NLOutgoingPsi(int outb_face_counter,
              int face_dof, int n)
{
  if (outb_face_counter>nonlocal_outb_face_deplocI_slot.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid number of outb_face_counter " << outb_face_counter
      << " max allowed " << nonlocal_outb_face_deplocI_slot.size();
    exit(EXIT_FAILURE);
  }

  int depLocI = nonlocal_outb_face_deplocI_slot[outb_face_counter].first;
  int slot    = nonlocal_outb_face_deplocI_slot[outb_face_counter].second;
  int nonlocal_psi_Gn_blockstride = deplocI_face_dof_count[depLocI];

  int index =
    nonlocal_psi_Gn_blockstride*G*n +
    slot*G + face_dof*G;

  if ((index<0) ||
      (index>ref_deplocI_outgoing_psi->operator[](depLocI).size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid index " << index
      << " encountered in non-local outgoing Psi"
      << " max allowed " << ref_deplocI_outgoing_psi->operator[](depLocI).size();
    exit(EXIT_FAILURE);
  }

  return &ref_deplocI_outgoing_psi->operator[](depLocI)[index];
}

//###################################################################
/**Given a sweep ordering index, the incoming face counter,
 * the incoming face dof, this function computes the location
 * where to store this position's outgoing psi and returns a reference
 * to it.*/
double*  chi_mesh::sweep_management::AUX_FLUDS::
UpwindPsi(int cell_so_index, int inc_face_counter,
          int face_dof,int g, int n)
{
  // Face category
  int fc = so_cell_inco_face_face_category[cell_so_index][inc_face_counter];

  if (fc >= 0)
  {
    size_t index =
      local_psi_Gn_block_strideG[fc]*n +
      so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].first*
      local_psi_stride[fc]*G +
      so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].
        second[face_dof]*G + g;

    return &(ref_local_psi->operator[](fc))[index];
  }
  else
  {
    size_t index =
      delayed_local_psi_Gn_block_strideG*n +
      so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].first*
      delayed_local_psi_stride*G +
      so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].
        second[face_dof]*G + g;

    return &ref_delayed_local_psi->operator[](index);
  }

}

//###################################################################
/**Given a sweep ordering index, the incoming face counter,
 * the incoming face dof, this function computes the location
 * where to obtain the position's upwind psi.*/
double*  chi_mesh::sweep_management::AUX_FLUDS::
NLUpwindPsi(int nonl_inc_face_counter,
            int face_dof,int g, int n)
{
  int prelocI =
    nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].first;

  if (prelocI>=0)
  {
    int nonlocal_psi_Gn_blockstride = prelocI_face_dof_count[prelocI];
    int slot =
      nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].second.first;

    int mapped_dof =
      nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].
        second.second[face_dof];

    int index =
      nonlocal_psi_Gn_blockstride*G*n +
      slot*G +
      mapped_dof*G + g;

    return &ref_prelocI_outgoing_psi->operator[](prelocI)[index];
  }
  else
  {
    prelocI =
      delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].first;

    int nonlocal_psi_Gn_blockstride = delayed_prelocI_face_dof_count[prelocI];
    int slot =
      delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].second.first;

    int mapped_dof =
      delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].
        second.second[face_dof];

    int index =
      nonlocal_psi_Gn_blockstride*G*n +
      slot*G +
      mapped_dof*G + g;

    return &ref_delayed_prelocI_outgoing_psi->operator[](prelocI)[index];
  }


}
