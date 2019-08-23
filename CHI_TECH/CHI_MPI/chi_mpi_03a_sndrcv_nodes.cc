#include "chi_mpi.h"

//###################################################################
/**This function sends a vector of nodes to all child processes*/
void CHI_MPI::BroadcastNodes(std::vector<chi_mesh::Vertex> *nodes)
{
  struct NODE_INFO
  {
    double xyz[3];
  };
  NODE_INFO* node_stack = new NODE_INFO[nodes->size()];
  for (unsigned v=0; v< nodes->size(); v++)
  {
    node_stack[v].xyz[0] = (*nodes)[v].x;
    node_stack[v].xyz[1] = (*nodes)[v].y;
    node_stack[v].xyz[2] = (*nodes)[v].z;
  }

  for (int k=1;k<this->process_count; k++)
  {
    MPI_Send(node_stack,nodes->size(),
             NODE_INFO_C, k,123,MPI_COMM_WORLD);
  }
  delete [] node_stack;
}



//###################################################################
/**This function receives nodes from the master.*/
void CHI_MPI::ReceiveNodes(std::vector<chi_mesh::Vertex> *nodes)
{
  //================================================== Find amount of nodes to
  //                                                   be received
  int node_count;
  MPI_Status status;
  MPI_Probe(0, 123, MPI_COMM_WORLD, &status);

  MPI_Get_count(&status, NODE_INFO_C, &node_count);

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
    chi_mesh::Node new_node;

    new_node.x = node_stack[k].xyz[0];
    new_node.y = node_stack[k].xyz[1];
    new_node.z = node_stack[k].xyz[2];

    nodes->push_back(new_node);
}
  delete [] node_stack;
}