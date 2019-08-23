#include "chi_mpi.h"
#include "../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../CHI_MESH/CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../CHI_MESH/CHI_CELL/cell.h"
#include "../CHI_MESH/CHI_CELL/cell_triangle.h"



//###################################################################
/**Broadcasts cells to child processes*/
void CHI_MPI::
  BroadcastCellSets()
{
  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Get master mesh continuum
  chi_mesh::CELL_SET* master_cell_set = mesh_handler->cell_sets[0];
  chi_mesh::MeshContinuum* mesh_continuum = master_cell_set->mesh_continuum;

  //================================================== Broadcast Nodes
  struct NODE_INFO
  {
    double xyz[3];
  };
  NODE_INFO* node_stack = new NODE_INFO[mesh_continuum->nodes.size()];
  for (unsigned v=0; v< mesh_continuum->nodes.size(); v++)
  {
   node_stack[v].xyz[0] = mesh_continuum->nodes[v]->x;
   node_stack[v].xyz[1] = mesh_continuum->nodes[v]->y;
   node_stack[v].xyz[2] = mesh_continuum->nodes[v]->z;
  }


  fprintf(stdout, "Broadcasting %d nodes\n",mesh_continuum->nodes.size());

  for (int k=1;k<this->process_count; k++)
  {
    MPI_Send(node_stack,mesh_continuum->nodes.size(),
             NODE_INFO_C, k,123,MPI_COMM_WORLD);
  }
  delete [] node_stack;

  //================================================== Broadcast item_id
  struct CELL_INFO
  {
    int v_index[3];
    int e_index0[4];
    int e_index1[4];
    int e_index2[4];
  };

  CELL_INFO* cell_stack = new CELL_INFO[mesh_continuum->cells.size()];
  for (unsigned c=0;c<mesh_continuum->cells.size();c++)
  {
    if (typeid(*mesh_continuum->cells[c]) == typeid(chi_mesh::CellTriangle))
    {
      cell_stack[c].v_index[0] = ((chi_mesh::CellTriangle*)mesh_continuum->
         cells[c])->v_index[0];
      cell_stack[c].v_index[1] = ((chi_mesh::CellTriangle*)mesh_continuum->
         cells[c])->v_index[1];
      cell_stack[c].v_index[2] = ((chi_mesh::CellTriangle*)mesh_continuum->
         cells[c])->v_index[2];

      cell_stack[c].e_index0[0] = ((chi_mesh::CellTriangle*)mesh_continuum->
       cells[c])->e_index[0][0];
      cell_stack[c].e_index0[1] = ((chi_mesh::CellTriangle*)mesh_continuum->
       cells[c])->e_index[0][1];
      cell_stack[c].e_index0[2] = ((chi_mesh::CellTriangle*)mesh_continuum->
       cells[c])->e_index[0][2];
      cell_stack[c].e_index0[3] = ((chi_mesh::CellTriangle*)mesh_continuum->
       cells[c])->e_index[0][3];

      cell_stack[c].e_index1[0] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[1][0];
      cell_stack[c].e_index1[1] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[1][1];
      cell_stack[c].e_index1[2] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[1][2];
      cell_stack[c].e_index1[3] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[1][3];

      cell_stack[c].e_index2[0] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[2][0];
      cell_stack[c].e_index2[1] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[2][1];
      cell_stack[c].e_index2[2] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[2][2];
      cell_stack[c].e_index2[3] = ((chi_mesh::CellTriangle*)mesh_continuum->
        cells[c])->e_index[2][3];
    }
    else
    {
      fprintf(stderr,"Wrong cell type in BroadcastCellSets\n");
      exit(EXIT_FAILURE);
    }
  }//for c

  fprintf(stdout, "Broadcasting %d item_id\n",mesh_continuum->cells.size());

  for (int k=1;k<this->process_count; k++)
  {
    MPI_Send(cell_stack,mesh_continuum->cells.size(),
             CELL_INFO_C, k,124,MPI_COMM_WORLD);
  }
  delete [] cell_stack;

  //================================================== Broadcast cell set
  //for (unsigned d=1; d<mesh_handler->cell_sets.size(); d++)
  for (int k=1;k<this->process_count; k++)
  {
    chi_mesh::CELL_SET* cell_set = mesh_handler->cell_sets[k];
    int* cell_ids = new int[cell_set->cells_allocated.size()];

    for (unsigned ci=0; ci< cell_set->cells_allocated.size(); ci++)
    {
      cell_ids[ci] = cell_set->cells_allocated[ci];
    }

    fprintf(stdout, "Broadcasting %d cell set id's\n",
            cell_set->cells_allocated.size());
    MPI_Send(cell_ids, cell_set->cells_allocated.size(),
             MPI_INT,k,125,MPI_COMM_WORLD);
  }

}