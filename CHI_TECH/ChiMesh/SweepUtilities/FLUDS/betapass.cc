#include "ChiMesh/SweepUtilities/chi_FLUDS.h"

#include "ChiMesh/SweepUtilities/chi_SPDS.h"

#include <ChiMesh/Cell/cell.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog     chi_log;
extern ChiMPI     chi_mpi;

//###################################################################
/**Receives and send predecessor data.*/
void chi_mesh::sweep_management::FLUDS::
InitializeBetaElements(chi_mesh::sweep_management::SPDS* spds, int tag_index)
{
  chi_mesh::MeshContinuum*         grid = spds->grid;
  chi_mesh::sweep_management::SPLS* spls = spds->spls;

//  chi_log.Log(LOG_0) << "Initializing FLUDS Beta elements";

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

    MPI_Recv(face_indices.data(),amount_to_receive,MPI_INT,
             locJ,101+tag_index,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    DeSerializeCellInfo(prelocI_cell_views[prelocI], &face_indices,
                        prelocI_face_dof_count[prelocI]);
  }

  //=============================================== Send successor information
  for (int deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
  {
    int locJ = spds->location_successors[deplocI];

    std::vector<int>::iterator delayed_successor =
    std::find(spds->delayed_location_successors.begin(),
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
  for (int csoi=0; csoi<spls->item_id.size(); csoi++)
  {
    int  cell_g_index = spls->item_id[csoi];
    auto cell         = grid->cells[cell_g_index];

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      TSlab* slab_cell = (TSlab*)cell;
      NonLocalIncidentMapping(slab_cell, spds);
    }//if slab
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      TPolygon* poly_cell = (TPolygon*)cell;
      NonLocalIncidentMapping(poly_cell, spds);
    }//if polyhedron
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      TPolyhedron* polyh_cell = (TPolyhedron*)cell;
      NonLocalIncidentMapping(polyh_cell, spds);
    }//if polyhedron

  }//for csoi

  //================================================== Clear unneccesary data
  auto empty_vector = std::vector<std::vector<CompactCellView>>(0);
  deplocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  prelocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  delayed_prelocI_cell_views.swap(empty_vector);
}