#include "FLUDS.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog&     chi_log;
extern ChiMPI&      chi_mpi;

#include <algorithm>

//###################################################################
/**Receives and send predecessor data.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
InitializeBetaElements(chi_mesh::sweep_management::SPDS* spds, int tag_index)
{
  chi_mesh::MeshContinuum*         grid = spds->grid;
  chi_mesh::sweep_management::SPLS& spls = spds->spls;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The first two major steps here are: Send delayed successor information
  // and Receive delayed predecessor information. The send portion is done
  // first because the delayed information does not follow the
  // Task Dependency Graph and hence when a location receives its delayed
  // information, the information might not have been sent yet.

  //=============================================== Send delayed successor information
  std::vector<MPI_Request>  send_requests;
  send_requests.resize(spds->location_successors.size(),MPI_Request());
  std::vector<std::vector<int>>
      multi_face_indices(spds->location_successors.size(),std::vector<int>());
  for (int deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
  {
    int locJ = spds->location_successors[deplocI];

    std::vector<int>::iterator delayed_successor =
        std::find(spds->delayed_location_successors.begin(),
                  spds->delayed_location_successors.end(),
                  locJ);
    if ((delayed_successor == spds->delayed_location_successors.end()))
      continue;

    std::vector<CompactCellView>* cell_views = &deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views,multi_face_indices[deplocI],
                      deplocI_face_dof_count[deplocI]);


    MPI_Isend(multi_face_indices[deplocI].data(),
              multi_face_indices[deplocI].size(),
              MPI_INT,locJ,101+tag_index,
              MPI_COMM_WORLD,&send_requests[deplocI]);

    //TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //=============================================== Receive delayed predecessor
  //                                                information
  delayed_prelocI_cell_views.resize(spds->delayed_location_dependencies.size(),
                                    std::vector<CompactCellView>());
  delayed_prelocI_face_dof_count.resize(spds->delayed_location_dependencies.size(),0);
  for (int prelocI=0; prelocI<spds->delayed_location_dependencies.size(); prelocI++)
  {
    int locJ = spds->delayed_location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ,101+tag_index,MPI_COMM_WORLD,&probe_status);

    int amount_to_receive=0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive );

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive,0);

    MPI_Recv(face_indices.data(),amount_to_receive,MPI_INT,
             locJ,101+tag_index,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    DeSerializeCellInfo(delayed_prelocI_cell_views[prelocI], &face_indices,
                        delayed_prelocI_face_dof_count[prelocI]);
  }

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The next two operations is to receive predecessor information followed
  // by the sending of successor information. The receives are blocking but
  // will cause a dead lock because the successor/predecessor combination
  // follows the TDG.

  //=============================================== Receive predecessor
  //                                                information
  prelocI_cell_views.resize(spds->location_dependencies.size(),
                            std::vector<CompactCellView>());
  prelocI_face_dof_count.resize(spds->location_dependencies.size(),0);
  for (int prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
  {
    int locJ = spds->location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ,101+tag_index,MPI_COMM_WORLD,&probe_status);

    int amount_to_receive=0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive );

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive,0);

    MPI_Recv(face_indices.data(),
             amount_to_receive,
             MPI_INT,
             locJ,101+tag_index,
             MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    DeSerializeCellInfo(prelocI_cell_views[prelocI], &face_indices,
                        prelocI_face_dof_count[prelocI]);
  }

  //=============================================== Send successor information
  for (int deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
  {
    int locJ = spds->location_successors[deplocI];

    auto delayed_successor = std::find(spds->delayed_location_successors.begin(),
                                       spds->delayed_location_successors.end(),
                                       locJ);
    if ((delayed_successor != spds->delayed_location_successors.end()))
      continue;

    std::vector<CompactCellView>* cell_views = &deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views,multi_face_indices[deplocI],
                      deplocI_face_dof_count[deplocI]);


    MPI_Isend(multi_face_indices[deplocI].data(),
              multi_face_indices[deplocI].size(),
              MPI_INT,locJ,101+tag_index,
              MPI_COMM_WORLD,&send_requests[deplocI]);

    //TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //================================================== Verify sends completed
  for (int deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
    MPI_Wait(&send_requests[deplocI],MPI_STATUS_IGNORE);
  multi_face_indices.clear();
  multi_face_indices.shrink_to_fit();


  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // In the next process we loop over cells in the sweep order and perform
  // the non-local face mappings. This is dependent on having the compact
  // cellviews on the partition interfaces.


  //================================================== Loop over cells in sorder
  for (int csoi=0; csoi<spls.item_id.size(); csoi++)
  {
    int cell_local_index = spls.item_id[csoi];
    auto cell = &grid->local_cells[cell_local_index];

    NonLocalIncidentMapping(cell, spds);
  }//for csoi

  deplocI_cell_views.clear();
  deplocI_cell_views.shrink_to_fit();

  prelocI_cell_views.clear();
  prelocI_cell_views.shrink_to_fit();

  delayed_prelocI_cell_views.clear();
  delayed_prelocI_cell_views.shrink_to_fit();

  //================================================== Clear unneccesary data
  auto empty_vector = std::vector<std::vector<CompactCellView>>(0);
  deplocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  prelocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  delayed_prelocI_cell_views.swap(empty_vector);
}

//###################################################################
/**This cell takes a hierarchy of a cell compact view and
 * serializes it for MPI transmission. This is easy since all
 * the values are integers.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
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

    size_t num_faces = cell_face_views->size();
    for (size_t f=0; f<num_faces; f++)
    {
      face_indices.push_back(glob_index);
      face_indices.push_back((*cell_face_views)[f].first);
      std::vector<int>* face_vertices = &(*cell_face_views)[f].second;

      size_t num_verts = face_vertices->size();
      for (int fi=0; fi<num_verts; fi++)
      {
        face_indices.push_back((*face_vertices)[fi]);
      }
    }
  }
}

//###################################################################
/**Deserializes face indices.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                    std::vector<int>* face_indices,
                    int& num_face_dofs)
{
  num_face_dofs = (*face_indices)[0];
  int num_cells     = (*face_indices)[1];

  cell_views.resize(num_cells);

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
        cell_views.emplace_back(); c++;
        cell_views[c].first = -entry-1;

        cell_views[c].second.emplace_back();f=0;

        v=0; last_cell = -entry;

        cell_views[c].second[f].first = (*face_indices)[k+1]; k++;
      } else
      {
        cell_views[c].second.emplace_back(); f++; v=0;

        cell_views[c].second[f].first = (*face_indices)[k+1]; k++;
      }
    }
      //================================= Face vertex
    else
    {
      cell_views[c].second[f].second.push_back(entry);
      v++;
    }
    k++;
  }//while k
}