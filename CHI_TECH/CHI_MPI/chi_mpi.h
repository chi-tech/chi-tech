#ifndef _chi_mpi_h
#define _chi_mpi_h

#include <mpi.h>
#include "../CHI_MESH/chi_mesh.h"

class CHI_MPI_COMMUNICATOR_SET
{
public:
  std::vector<MPI_Comm>  communicators;
  std::vector<MPI_Group> location_groups;
  MPI_Group              world_group;

public:
  int MapIonJ(int locI, int locJ)
  {
    int group_rank;
    MPI_Group_translate_ranks(world_group,1,&locI,
                              location_groups[locJ],&group_rank);

    return group_rank;
  }
};

//################################################################### Class def
/**An object for storing various MPI states.*/
class CHI_MPI
{
public:
  int location_id;
  int process_count;
  MPI_Datatype NODE_INFO_C;
  MPI_Datatype TRIFACE_INFO_C;
  MPI_Datatype CELL_INFO_C;

  MPI_Datatype LOC_SWP_DEP_C;

public:
  CHI_MPI()
  {
    location_id = 0;
    process_count = 1;
  }
  //01
  void Initialize();

  //02
  void BroadcastCellSets();
  void ReceiveCellSets();

  //03a
  void BroadcastNodes(std::vector<chi_mesh::Vertex>* nodes);
  void ReceiveNodes(std::vector<chi_mesh::Vertex>* nodes);

  //03a
  void BroadcastTriFaces(std::vector<chi_mesh::Face>* faces);
  void ReceiveTriFaces(std::vector<chi_mesh::Face>* faces);

  //03c
  void SendSweepDependency(int dest, std::vector<int>* dependencies);
  void ReceiveSweepDependency(int sorc, std::vector<int>* dependencies);
};




#endif