#ifndef CHI_MPI_H
#define CHI_MPI_H

#include <mpi.h>
#include "../ChiMesh/chi_mesh.h"

//################################################################### Class def
/**Simple implementation a communicator set.*/
class ChiMPICommunicatorSet
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
class ChiMPI
{
private:
  int m_location_id = 0;
  int m_process_count = 1;

  bool m_location_id_set = false;
  bool m_process_count_set = false;

public:
  const int& location_id = m_location_id;
  const int& process_count = m_process_count;

private:
  static ChiMPI instance;
  ChiMPI() = default;

public:
  void SetLocationID(int in_location_id)
  {
    if (not m_location_id_set)
      m_location_id = in_location_id;
    m_location_id_set = true;
  }
  void SetProcessCount(int in_process_count)
  {
    if (not m_process_count_set)
      m_process_count = in_process_count;
    m_process_count_set = true;
  }
  static ChiMPI& GetInstance() noexcept {return instance;}
};




#endif