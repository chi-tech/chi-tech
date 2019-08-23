#include "chi_mpi.h"

//###################################################################
/**This function sends a vector of dependency indices to the
 * specified location*/
void CHI_MPI::SendSweepDependency(int dest, std::vector<int> *dependencies)
{
  struct DEPENDECY_INFO
  {
    int deps[50];
    int length;
    int location;
  };

  DEPENDECY_INFO dep_info;
  dep_info.location = this->location_id;
  dep_info.length = dependencies->size();
  for (int k=0; k<dep_info.length; k++)
  {
    if (k>49)
    {
      fprintf(stderr,"ERROR: SendSweepDependency. Number of dependcies "
                     "greater than allowed.\n");
      exit(EXIT_FAILURE);
    }
    dep_info.deps[k] = (*dependencies)[k];
  }

  MPI_Send(&dep_info,1,LOC_SWP_DEP_C,dest,124,MPI_COMM_WORLD);
}

//###################################################################
/**This function receives a vector of dependcy indices from the
 * specified location.*/
void CHI_MPI::ReceiveSweepDependency(int sorc, std::vector<int> *dependencies)
{
  //============================================= Check for message
  int message_count;
  MPI_Status status;
  MPI_Probe(sorc,124,MPI_COMM_WORLD,&status);

  MPI_Get_count(&status, LOC_SWP_DEP_C, &message_count);

  //============================================= Receive the information
  struct DEPENDECY_INFO
  {
    int deps[50];
    int length;
    int location;
  };

  DEPENDECY_INFO dep_info;

  MPI_Recv(&dep_info,1,LOC_SWP_DEP_C,sorc,124,MPI_COMM_WORLD,&status);

  for (int k=0; k<dep_info.length; k++)
  {
    if (k>49)
    {
      fprintf(stderr,"ERROR: RecvSweepDependency. Number of dependcies "
                     "greater than allowed.\n");
      exit(EXIT_FAILURE);
    }
    dependencies->push_back(dep_info.deps[k]);
  }
}