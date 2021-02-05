#ifndef _chi_spds_h
#define _chi_spds_h

#include "ChiMesh/SweepUtilities/SPLS/SPLS.h"

#include <memory>

namespace chi_mesh::sweep_management
{
  struct SPDS;
//  typedef std::shared_ptr<SPDS> SPDS_ptr;
}

//###################################################################
/**Contains multiple levels*/
struct chi_mesh::sweep_management::SPDS
{
  double                   polar;
  double                   azimuthal;
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

  //======================================== Default constructor
  SPDS()
  {  }

  int MapLocJToPrelocI(int locJ);
  int MapLocJToDeplocI(int locJ);

  void BuildTaskDependencyGraph(bool cycle_allowance_flag);
};

#endif