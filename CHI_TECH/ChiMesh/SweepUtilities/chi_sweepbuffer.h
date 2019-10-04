#ifndef _chi_sweepbuffer_h
#define _chi_sweepbuffer_h

#include "chi_sweep.h"
#include <chi_mpi.h>

#ifndef FLAG_FINISHED
  #define FLAG_FINISHED     true
  #define FLAG_NOT_FINISHED false
#endif

typedef unsigned long long int u_ll_int;

namespace chi_mesh::SweepManagement
{

struct BoundaryTypes
{
  static const int VACUUM = 1;
  static const int INCIDENT_ISOTROPIC = 2;
};


//###################################################################
/**Handles the swift communication of interprocess communication
 * related to sweeping.*/
class SweepBuffer
{
private:
  chi_mesh::SweepManagement::AngleSet* angleset;
  bool initialized;
  bool done_sending;
  bool data_initialized;
  bool upstream_data_initialized;

  u_ll_int EAGER_LIMIT = 32000;

  std::vector<int> prelocI_message_count;
  std::vector<int> deplocI_message_count;
  std::vector<int> delayed_prelocI_message_count;

  std::vector<std::vector<u_ll_int>> prelocI_message_size;
  std::vector<std::vector<u_ll_int>> deplocI_message_size;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_size;

  std::vector<std::vector<u_ll_int>> prelocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> deplocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_blockpos;

  std::vector<std::vector<bool>> prelocI_message_available;
  std::vector<std::vector<bool>> deplocI_message_sent; //Might be redundant

  std::vector<std::vector<bool>> delayed_prelocI_message_available;


  std::vector<std::vector<MPI_Request>> deplocI_message_request;

  ChiMPICommunicatorSet*   comm_set;


public:
  SweepBuffer(chi_mesh::SweepManagement::AngleSet* ref_angleset,
              int sweep_eager_limit,
              ChiMPICommunicatorSet* in_comm_set);
  void CheckInitialized();
  void InitializeBuffers();
  void SendDownstreamPsi(int angle_set_num);
  void ReceiveDelayedData(int angle_set_num);
  void CheckDownstreamBuffersClear();
  bool CheckUpstreamPsiAvailable(int angle_set_num);
  void ClearReceiveBuffers();

  void Reset()
  {
    done_sending = false;
    data_initialized = false;
    upstream_data_initialized = false;

    for (int prelocI=0; prelocI<prelocI_message_available.size(); prelocI++)
    {
      for (int m=0; m<prelocI_message_available[prelocI].size(); m++)
      {
        prelocI_message_available[prelocI][m] = false;
      }
    }
  }

};
}
#endif