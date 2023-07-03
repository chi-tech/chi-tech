#ifndef CHI_SWEEPBUFFER_H
#define CHI_SWEEPBUFFER_H

#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "chi_mpi.h"

typedef unsigned long long int u_ll_int;

namespace chi
{
  class ChiMPICommunicatorSet;
}

namespace chi_mesh::sweep_management
{

//###################################################################
/**Handles the swift communication of interprocess communication
 * related to sweeping.*/
class SweepBuffer
{
private:
  chi_mesh::sweep_management::AngleSet* const angleset;
  const chi::ChiMPICommunicatorSet&   comm_set;

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

  std::vector<std::vector<bool>> prelocI_message_received;
  std::vector<std::vector<bool>> delayed_prelocI_message_received;

  std::vector<std::vector<MPI_Request>> deplocI_message_request;



public:
  int max_num_mess;

  SweepBuffer(chi_mesh::sweep_management::AngleSet* ref_angleset,
              int sweep_eager_limit,
              const chi::ChiMPICommunicatorSet& in_comm_set);
  bool DoneSending() const;
  void BuildMessageStructure();
  void InitializeDelayedUpstreamData();
  void InitializeLocalAndDownstreamBuffers();
  void SendDownstreamPsi(int angle_set_num);
  bool ReceiveDelayedData(int angle_set_num);
  void ClearDownstreamBuffers();
  AngleSetStatus ReceiveUpstreamPsi(int angle_set_num);
  void ClearLocalAndReceiveBuffers();
  void Reset();

};
}
#endif //CHI_SWEEPBUFFER_H