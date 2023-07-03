#ifndef CHI_SPDS_H
#define CHI_SPDS_H

#include "ChiMesh/SweepUtilities/SPLS/SPLS.h"

#include <memory>

namespace chi_mesh::sweep_management
{
  struct SPDS;
}

//###################################################################
/**Contains multiple levels*/
struct chi_mesh::sweep_management::SPDS
{
  chi_mesh::Vector3        omega;

  chi_mesh::MeshContinuumPtr grid;

  SPLS                     spls;
  std::vector<STDG>        global_sweep_planes;  ///< Processor sweep planes
  std::vector<int>         location_dependencies;
  std::vector<int>         location_successors;
  std::vector<int>         delayed_location_dependencies;
  std::vector<int>         delayed_location_successors;

  std::vector<std::pair<int,int>> local_cyclic_dependencies;

  std::vector<std::vector<int>> global_dependencies;

  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;

  //======================================== Default constructor
  SPDS() = default;

  int MapLocJToPrelocI(int locJ) const;
  int MapLocJToDeplocI(int locJ) const;

  void BuildTaskDependencyGraph(bool cycle_allowance_flag);
};

#endif //CHI_SPDS_H