#include "chi_mpi.h"
#include "../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../CHI_MESH/CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../CHI_MESH/CHI_CELL/cell_triangle.h"



//###################################################################
/**Broadcasts cells to child processes*/
void CHI_MPI::
ReceiveCellSets()
{
  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Create a cell set
  //                                                   into which receive
  chi_mesh::CELL_SET* cell_set = new chi_mesh::CELL_SET;
  cell_set->mesh_continuum = new chi_mesh::MeshContinuum;


  //================================================== Find amount of nodes to
  //                                                   be received
  int node_count;
  MPI_Status status;
  MPI_Probe(0, 123, MPI_COMM_WORLD, &status);

  MPI_Get_count(&status, NODE_INFO_C, &node_count);

  fprintf(stdout,"Process %d is receiving %d nodes\n",
                 this->location_id,node_count);

  //================================================== Receive nodes
  struct NODE_INFO
  {
    double xyz[3];
  };

  NODE_INFO* node_stack = new NODE_INFO[node_count];

  MPI_Recv(node_stack,node_count,
           NODE_INFO_C, 0,123,MPI_COMM_WORLD,&status);

  for (int k=0;k<node_count;k++)
  {
    chi_mesh::Node* new_node = new chi_mesh::Node;

    new_node->x = node_stack[k].xyz[0];
    new_node->x = node_stack[k].xyz[1];
    new_node->x = node_stack[k].xyz[2];

    cell_set->mesh_continuum->nodes.push_back(new_node);
  }

  //================================================== Find amount of item_id to
  //                                                   be received
  int cell_count;
  MPI_Probe(0, 124, MPI_COMM_WORLD, &status);

  MPI_Get_count(&status, CELL_INFO_C, &cell_count);

  fprintf(stdout,"Process %d is receiving %d item_id\n",
          this->location_id,cell_count);

  //================================================== Receive item_id
  struct CELL_INFO
  {
    int v_index[3];
    int e_index0[4];
    int e_index1[4];
    int e_index2[4];
  };

  CELL_INFO* cell_stack = new CELL_INFO[cell_count];

  MPI_Recv(cell_stack,cell_count,
           CELL_INFO_C, 0,124,MPI_COMM_WORLD,&status);

  for (int k=0;k<cell_count;k++)
  {
    //fprintf(stdout,"Cell %d\n",k);
    chi_mesh::CellTriangle* new_cell = new chi_mesh::CellTriangle;

    CELL_INFO* cell_info = &cell_stack[k];
    for (int e=0;e<3;e++)
    {
      new_cell->v_index[e] = cell_info->v_index[e];
    }

    for (int e=0;e<4;e++)
    {
      new_cell->e_index[0][e] = cell_info->e_index0[e];
    }

    for (int e=0;e<4;e++)
    {
      new_cell->e_index[1][e] = cell_info->e_index1[e];
    }

    for (int e=0;e<4;e++)
    {
      new_cell->e_index[2][e] = cell_info->e_index2[e];
    }

    cell_set->mesh_continuum->cells.push_back(new_cell);
  }

  //================================================== Receive cell set
  int cellid_count;
  MPI_Probe(0, 125, MPI_COMM_WORLD, &status);

  MPI_Get_count(&status, MPI_INT, &cellid_count);

  fprintf(stdout,"Process %d is receiving %d cell id's\n",
          this->location_id,cellid_count);

  int* cell_ids = new int[cellid_count];


  MPI_Recv(cell_ids,cellid_count,
           MPI_INT, 0,125,MPI_COMM_WORLD,&status);


  //============================================= Parsing the information
  for (unsigned ci=0; ci< cellid_count; ci++)
  {
    cell_set->cells_allocated.push_back(cell_ids[ci]);
  }

  mesh_handler->cell_sets.push_back(cell_set);




}