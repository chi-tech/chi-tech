#ifndef _chi_angleset_h
#define _chi_angleset_h

#include "../chi_mesh.h"
#include "chi_sweepbuffer.h"
#include "chi_sweep_boundaries.h"
#include <ChiMesh/SweepUtilities/chi_FLUDS.h>

#include <chi_mpi.h>

typedef chi_mesh::SweepManagement::BoundaryBase SweepBndry;

//###################################################################
/**Manages the workstages of a single angle set.*/
class chi_mesh::SweepManagement::AngleSet
{
private:
  int              num_grps;
  SPDS*            spds;
  bool             executed;

  chi_mesh::SweepManagement::SweepBuffer sweep_buffer;

public:
  FLUDS*                                 fluds;
  std::vector<int>                       angles;
  std::vector<SweepBndry*>&              ref_boundaries;
  int                                    ref_subset;

  //FLUDS
  std::vector<double>               local_psi;
  std::vector<double>               delayed_local_psi;
  std::vector<std::vector<double>>  deplocI_outgoing_psi;
  std::vector<std::vector<double>>  prelocI_outgoing_psi;
  std::vector<std::vector<double>>  boundryI_incoming_psi;

  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi;
  std::vector<double>               delayed_prelocI_norm;

  AngleSet(int in_numgrps,
           int in_ref_subset,
           SPDS* in_spds,
           std::vector<int>& angle_indices,
           std::vector<SweepBndry*>& sim_boundaries,
           int sweep_eager_limit,
           ChiMPICommunicatorSet* in_comm_set):
           sweep_buffer(this,sweep_eager_limit,in_comm_set),
           ref_boundaries(sim_boundaries)
  {
    num_grps = in_numgrps;
    spds     = in_spds;
    executed = false;
    ref_subset = in_ref_subset;
    std::copy(angle_indices.begin(),
              angle_indices.end(),
              std::back_inserter(angles));

    fluds = new chi_mesh::SweepManagement::FLUDS(num_grps);
    fluds->InitializeAlphaElements(spds);
    fluds->InitializeBetaElements(spds);
  };

  void InitializeDelayedUpstreamData();

  SPDS* GetSPDS()
  {
    return spds;
  }

  int GetNumGrps() {return num_grps;}

  void EnsureClearedBuffers() {
    if (executed) sweep_buffer.CheckDownstreamBuffersClear();}




  bool AngleSetAdvance(chi_mesh::SweepManagement::SweepChunk *sweep_chunk,
                       int angle_set_num);
  void ResetSweepBuffers()
  {
    sweep_buffer.Reset();
    executed = false;
  }
  void ReceiveDelayedData(int angle_set_num)
  {
    sweep_buffer.ReceiveDelayedData(angle_set_num);
  }

  double* PsiBndry(int bndry_face_count, int bndry_map,
                      int face_dof, int g,int angle_num)
  {
    double* Psi = &ref_boundaries[bndry_map]->boundary_flux.data()[g];
    return Psi;
  }

};

#endif