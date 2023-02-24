#ifndef CHITECH_CHI_MPI_COMMSET_H
#define CHITECH_CHI_MPI_COMMSET_H

#include <mpi.h>
#include "../ChiMesh/chi_mesh.h"
#include "chi_runtime.h"

namespace chi_objects
{

//################################################################### Class def
/**Simple implementation a communicator set.*/
class ChiMPICommunicatorSet
{
private:
  std::vector<MPI_Comm>  communicators_;
  std::vector<MPI_Group> location_groups_;
  MPI_Group              world_group_;

public:
  ChiMPICommunicatorSet(std::vector<MPI_Comm>&  communicators,
                        std::vector<MPI_Group>& location_groups,
                        MPI_Group&               world_group) :
    communicators_(communicators),
    location_groups_(location_groups),
    world_group_(world_group)
  {}

  MPI_Comm LocICommunicator(int locI) const
  {
    return communicators_[locI];
  }

  int MapIonJ(int locI, int locJ) const
  {
    int group_rank;
    MPI_Group_translate_ranks(world_group_, 1, &locI,
                              location_groups_[locJ], &group_rank);

    return group_rank;
  }
};
}//namespace chi_objects




#endif //CHITECH_CHI_MPI_COMMSET_H
