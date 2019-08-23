#include "chi_mpi.h"

//###################################################################
/**This function broadcast a set of triangle faces to all child
 * processes.*/
void CHI_MPI::BroadcastTriFaces(std::vector<chi_mesh::Face> *faces)
{
  struct FACE_INFO
  {
    int v_index[3];
    int e_index0[4];
    int e_index1[4];
    int e_index2[4];
  };

  FACE_INFO* face_stack = new FACE_INFO[faces->size()];
  for (unsigned c=0;c<faces->size();c++)
  {
    face_stack[c].v_index[0] = (*faces)[c].v_index[0];
    face_stack[c].v_index[1] = (*faces)[c].v_index[1];
    face_stack[c].v_index[2] = (*faces)[c].v_index[2];

    face_stack[c].e_index0[0] = (*faces)[c].e_index[0][0];
    face_stack[c].e_index0[1] = (*faces)[c].e_index[0][1];
    face_stack[c].e_index0[2] = (*faces)[c].e_index[0][2];
    face_stack[c].e_index0[3] = (*faces)[c].e_index[0][3];

    face_stack[c].e_index1[0] = (*faces)[c].e_index[1][0];
    face_stack[c].e_index1[1] = (*faces)[c].e_index[1][1];
    face_stack[c].e_index1[2] = (*faces)[c].e_index[1][2];
    face_stack[c].e_index1[3] = (*faces)[c].e_index[1][3];

    face_stack[c].e_index2[0] = (*faces)[c].e_index[2][0];
    face_stack[c].e_index2[1] = (*faces)[c].e_index[2][1];
    face_stack[c].e_index2[2] = (*faces)[c].e_index[2][2];
    face_stack[c].e_index2[3] = (*faces)[c].e_index[2][3];
  }//for c

  for (int k=1;k<this->process_count; k++)
  {
    MPI_Send(face_stack,faces->size(),
             TRIFACE_INFO_C, k,124,MPI_COMM_WORLD);
  }
  delete [] face_stack;

}


//###################################################################
/**This function receives a set of triangle faces from the master.*/
void CHI_MPI::ReceiveTriFaces(std::vector<chi_mesh::Face> *faces)
{
//================================================== Find amount of item_id to
  //                                                   be received
  int face_count;
  MPI_Status status;
  MPI_Probe(0, 124, MPI_COMM_WORLD, &status);

  MPI_Get_count(&status, TRIFACE_INFO_C, &face_count);

  //================================================== Receive item_id
  struct FACE_INFO
  {
    int v_index[3];
    int e_index0[4];
    int e_index1[4];
    int e_index2[4];
  };

  FACE_INFO* face_stack = new FACE_INFO[face_count];

  MPI_Recv(face_stack,face_count,
           TRIFACE_INFO_C, 0,124,MPI_COMM_WORLD,&status);

  for (int k=0;k<face_count;k++)
  {
    chi_mesh::Face new_face;

    FACE_INFO* face_info = &face_stack[k];

    for (int e=0;e<3;e++)
    {
      new_face.v_index[e] = face_info->v_index[e];
    }

    for (int e=0;e<4;e++)
    {
      new_face.e_index[0][e] = face_info->e_index0[e];
    }

    for (int e=0;e<4;e++)
    {
      new_face.e_index[1][e] = face_info->e_index1[e];
    }

    for (int e=0;e<4;e++)
    {
      new_face.e_index[2][e] = face_info->e_index2[e];
    }

    faces->push_back(new_face);
  }

}