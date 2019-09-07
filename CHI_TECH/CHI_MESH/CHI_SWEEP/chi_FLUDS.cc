#include "chi_FLUDS.h"

#include "../../CHI_CONSOLE/chi_console.h"
#include "chi_SPDS.h"

#include "../../CHI_MESH/CHI_CELL/cell.h"
#include "../../CHI_MESH/CHI_CELL/cell_slab.h"
#include "../../CHI_MESH/CHI_CELL/cell_polygon.h"
#include "../../CHI_MESH/CHI_CELL/cell_polyhedron.h"

#include <chi_log.h>
#include <chi_mpi.h>

#include <iomanip>

extern CHI_CONSOLE chi_console;
extern ChiLog     chi_log;
extern ChiMPI     chi_mpi;

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector
 * and returns a reference to it.*/
double*  chi_mesh::SweepManagement::FLUDS::
 OutgoingPsi(int cell_so_index, int outb_face_counter,
              int face_dof, int n)
{
  int index =
    local_psi_Gn_block_stride*G*n +
    so_cell_outb_face_slot_indices[cell_so_index][outb_face_counter]*
    local_psi_stride*G +
    face_dof*G;

  return &ref_local_psi->operator[](index);
}

//###################################################################
/**Given a */
double*  chi_mesh::SweepManagement::FLUDS::
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
double*  chi_mesh::SweepManagement::FLUDS::
UpwindPsi(int cell_so_index, int inc_face_counter,
             int face_dof,int g, int n)
{
  int index =
    local_psi_Gn_block_stride*G*n +
    so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].first*
    local_psi_stride*G +
    so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].
      second[face_dof]*G + g;

  return &ref_local_psi->operator[](index);
}

//###################################################################
/**Given a sweep ordering index, the incoming face counter,
 * the incoming face dof, this function computes the location
 * where to obtain the position's upwind psi.*/
double*  chi_mesh::SweepManagement::FLUDS::
 NLUpwindPsi(int nonl_inc_face_counter,
          int face_dof,int g, int n)
{
  int prelocI =
    nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].first;

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

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector.*/
void  chi_mesh::SweepManagement::FLUDS::
AddFaceViewToDepLocI(int deplocI, int cell_g_index, int face_slot,
                     TVertexFace face_v_index)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (int c=0; c<deplocI_cell_views[deplocI].size(); c++)
  {
    if (deplocI_cell_views[deplocI][c].first == cell_g_index)
    {
      cell_already_there = true;
      deplocI_cell_views[deplocI][c].second.
        push_back(  CompactFaceView(face_slot,std::vector<int>(1,face_v_index)));
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.
      push_back(  CompactFaceView(face_slot,std::vector<int>(1,face_v_index)));

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }
}

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector.*/
void  chi_mesh::SweepManagement::FLUDS::
AddFaceViewToDepLocI(int deplocI, int cell_g_index, int face_slot,
                     TEdgeFace edge_v_indices)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (int c=0; c<deplocI_cell_views[deplocI].size(); c++)
  {
    if (deplocI_cell_views[deplocI][c].first == cell_g_index)
    {
      cell_already_there = true;
      std::vector<int> verts;
      verts.push_back(edge_v_indices[0]);
      verts.push_back(edge_v_indices[1]);


      deplocI_cell_views[deplocI][c].second.
        push_back(  CompactFaceView(face_slot,verts));
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    std::vector<int> verts;
    verts.push_back(edge_v_indices[0]);
    verts.push_back(edge_v_indices[1]);

    new_cell_view.second.
      push_back(  CompactFaceView(face_slot,verts));

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }
}

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector.*/
void  chi_mesh::SweepManagement::FLUDS::
  AddFaceViewToDepLocI(int deplocI, int cell_g_index, int face_slot,
                       TPolyFace *poly_face)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (int c=0; c<deplocI_cell_views[deplocI].size(); c++)
  {
    if (deplocI_cell_views[deplocI][c].first == cell_g_index)
    {
      cell_already_there = true;
      deplocI_cell_views[deplocI][c].second.
        push_back(  CompactFaceView(face_slot,poly_face->v_indices)  );
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.
      push_back(  CompactFaceView(face_slot,poly_face->v_indices)  );

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }


}








//###################################################################
/**This cell takes a hierarchy of a cell compact view and
 * serializes it for MPI transmission. This is easy since all
 * the values are integers.*/
void chi_mesh::SweepManagement::FLUDS::
SerializeCellInfo(std::vector<CompactCellView>* cell_views,
                  std::vector<int>& face_indices,
                  int num_face_dofs)
{
  int num_cells = cell_views->size();

  //======================== First entry is number of face dofs
  face_indices.push_back(num_face_dofs);

  //======================== Second entry is amount of cells
  face_indices.push_back(num_cells);

  //======================== Third entry is negative global cell index
  // Each time a negative entry occurs it denotes a cell face but
  // the actual number is -cell_g_index-1. The offset is necessary
  // for evaluating the negative. The offset is restored during the
  // deserialization process.
  // It is followed by a positive number which is the store location
  // of the face
  for (int c=0; c<num_cells; c++)
  {
    int glob_index = -(*cell_views)[c].first-1;

    std::vector<CompactFaceView>* cell_face_views =
      &(*cell_views)[c].second;

    int num_faces = cell_face_views->size();
    for (int f=0; f<num_faces; f++)
    {
      face_indices.push_back(glob_index);
      face_indices.push_back((*cell_face_views)[f].first);
      std::vector<int>* face_vertices = &(*cell_face_views)[f].second;

      int num_verts = face_vertices->size();
      for (int fi=0; fi<num_verts; fi++)
      {
        face_indices.push_back((*face_vertices)[fi]);
      }
    }
  }
}

//###################################################################
/**Deserializes face indices.*/
void chi_mesh::SweepManagement::FLUDS::
DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                    std::vector<int>* face_indices,
                    int& num_face_dofs)
{
  num_face_dofs = (*face_indices)[0];
  int num_cells     = (*face_indices)[1];

  cell_views.resize(num_cells);
  //chi_log.Log(LOG_ALL) << "Number of cells= " << num_cells;

  int k         =  2;
  int last_cell = -1;
  int c         = -1; //cell counter
  int f         = -1;
  int v         = -1;
  while (k<face_indices->size())
  {
    int entry = (*face_indices)[k];
    //================================= Cell/Face indicator
    if (entry < 0)
    {
      if (-entry != last_cell)
      {
        cell_views.push_back(CompactCellView()); c++;
        cell_views[c].first = -entry-1;

        cell_views[c].second.push_back(CompactFaceView());f=0;

        v=0; last_cell = -entry;

        cell_views[c].second[f].first = (*face_indices)[k+1]; k++;
      } else
      {
        cell_views[c].second.push_back(CompactFaceView()); f++; v=0;

        cell_views[c].second[f].first = (*face_indices)[k+1]; k++;
      }
    }
      //================================= Face vertex
    else
    {
      cell_views[c].second[f].second.push_back(entry);


//      chi_log.Log(LOG_ALL) << "Cell " << c
//                                << "(" << cell_views[c].first << ")"
//                                << " Face " << f
//                                << " Store " << cell_views[c].second[f].first
//                                << " Vert " << v
//                                << " Val= " << entry;


      v++;
    }
    k++;
  }//while k
}












