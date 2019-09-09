#ifndef _chi_angleset_h
#define _chi_angleset_h

#include "../chi_mesh.h"
#include "chi_sweepbuffer.h"

#include <chi_mpi.h>

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
  std::vector<std::pair<int,int>>*       ref_boundary_types;
  std::vector<std::vector<double>>*      ref_incident_P0_mg_boundaries;
  int                                    ref_subset;

  AngleSet(int in_numgrps,
           int in_ref_subset,
           SPDS* in_spds,
           std::vector<std::pair<int,int>>* sim_boundary_types,
           std::vector<std::vector<double>>* incident_P0_mg_boundaries,
           int sweep_eager_limit,
           ChiMPICommunicatorSet* in_comm_set):
           sweep_buffer(this,sweep_eager_limit,in_comm_set)
  {
    num_grps = in_numgrps;
    spds     = in_spds;
    ref_boundary_types = sim_boundary_types;
    ref_incident_P0_mg_boundaries = incident_P0_mg_boundaries;
    executed = false;
    ref_subset = in_ref_subset;
  };

  SPDS* GetSPDS()
  {
    return spds;
  }

  int GetNumGrps() {return num_grps;}

  void EnsureClearedBuffers() {
    if (executed) sweep_buffer.CheckDownstreamBuffersClear();}

  //FLUDS
  std::vector<double>               local_psi;
  std::vector<std::vector<double>>  deplocI_outgoing_psi;
  std::vector<std::vector<double>>  prelocI_outgoing_psi;
  std::vector<std::vector<double>>  boundryI_incoming_psi;


  bool AngleSetAdvance(chi_mesh::SweepManagement::SweepChunk *sweep_chunk,
                       int angle_set_num);
  void ResetSweepBuffers()
  {
    sweep_buffer.Reset();
    executed = false;
  }

  double* PsiBndry(int bndry_face_count, int bndry_map,
                      int face_dof, int g,int angle_num)
  {
    double* Psi = &boundryI_incoming_psi[bndry_map].data()[g];
    return Psi;
  }

};

#endif