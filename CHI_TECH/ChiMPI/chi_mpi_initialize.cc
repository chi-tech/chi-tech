#include "chi_mpi.h"
#include "stddef.h"

//###################################################################
/** This function initializes MPI datastructures.*/
void ChiMPI::Initialize()
{
  MPI_Datatype* data_types;
  int*          block_lengths;
  MPI_Aint*     block_displacements;


  //============================================= NODE INFO
  data_types          = new MPI_Datatype[1];
  block_lengths       = new int[1];
  block_displacements = new MPI_Aint[1];

  data_types[0]          = MPI_DOUBLE;   //Array type
  block_lengths[0]       = 3;            //Array size x,y,z
  block_displacements[0] = 0;            //Displacement between arrays

  MPI_Type_create_struct(1,block_lengths,
                         block_displacements,
                         data_types,
                         &NODE_INFO_C);
  MPI_Type_commit(&NODE_INFO_C);

  delete [] data_types;
  delete [] block_lengths;
  delete [] block_displacements;

  //============================================= TRIFACE INFO
  struct FACE_INFO
  {
    int v_index[3];
    int e_index0[4];
    int e_index1[4];
    int e_index2[4];
  };
  data_types          = new MPI_Datatype[4];
  block_lengths       = new int[4];
  block_displacements = new MPI_Aint[4];

  data_types[0]          = MPI_INT;
  data_types[1]          = MPI_INT;
  data_types[2]          = MPI_INT;
  data_types[3]          = MPI_INT;
  block_lengths[0]       = 3;
  block_lengths[1]       = 4;
  block_lengths[2]       = 4;
  block_lengths[3]       = 4;
  block_displacements[0] = offsetof(FACE_INFO,v_index);
  block_displacements[1] = offsetof(FACE_INFO,e_index0);
  block_displacements[2] = offsetof(FACE_INFO,e_index1);
  block_displacements[3] = offsetof(FACE_INFO,e_index2);

  MPI_Type_create_struct(4,block_lengths,
                         block_displacements,
                         data_types,
                         &TRIFACE_INFO_C);
  MPI_Type_commit(&TRIFACE_INFO_C);

  delete [] data_types;
  delete [] block_lengths;
  delete [] block_displacements;

  //============================================= CELL INFO
  data_types          = new MPI_Datatype[4];
  block_lengths       = new int[4];
  block_displacements = new MPI_Aint[4];

  data_types[0]          = MPI_INT;
  data_types[1]          = MPI_INT;
  data_types[2]          = MPI_INT;
  data_types[3]          = MPI_INT;
  block_lengths[0]       = 3;
  block_lengths[1]       = 4;
  block_lengths[2]       = 4;
  block_lengths[3]       = 4;
  block_displacements[0] = 0;
  block_displacements[1] = 3*sizeof(int);
  block_displacements[2] = 7*sizeof(int);
  block_displacements[3] = 11*sizeof(int);

  MPI_Type_create_struct(4,block_lengths,
                         block_displacements,
                         data_types,
                         &CELL_INFO_C);
  MPI_Type_commit(&CELL_INFO_C);

  delete [] data_types;
  delete [] block_lengths;
  delete [] block_displacements;

  //============================================= LOCATION SWEEP DEPENDENCY
  //We reserve a message size upto 50 sweep dependencies
  //for a given location. This is overkill for block
  //partitioning which could at max have 6 dependencies
  //but I figured I could support more for future
  //compatibility.
  //We also pass an integer indicating how much of
  //of this array of values are useful.
  data_types          = new MPI_Datatype[2];
  block_lengths       = new int[2];
  block_displacements = new MPI_Aint[2];

  data_types[0]          = MPI_INT;
  data_types[1]          = MPI_INT;
  block_lengths[0]       = 50;
  block_lengths[1]       = 2;
  block_displacements[0] = 0;
  block_displacements[1] = 50*sizeof(int);

  MPI_Type_create_struct(2,block_lengths,
                         block_displacements,
                         data_types,
                         &LOC_SWP_DEP_C);
  MPI_Type_commit(&LOC_SWP_DEP_C);

  delete [] data_types;
  delete [] block_lengths;
  delete [] block_displacements;

}