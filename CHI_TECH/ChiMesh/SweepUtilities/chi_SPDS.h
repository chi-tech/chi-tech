#ifndef _chi_spds_h
#define _chi_spds_h

#include "chi_SPLS.h"

//###################################################################
/**Contains multiple levels*/
struct chi_mesh::SweepManagement::SPDS
{
  double                   polar;
  double                   azimuthal;
  chi_mesh::Vector         omega;

  chi_mesh::MeshContinuum* grid;

  SPLS*                    spls;
  std::vector<STDG*>       global_sweep_planes;  ///< Processor sweep planes
  std::vector<int>         location_dependencies;
  std::vector<int>         location_successors;
  std::vector<int>         delayed_location_dependencies;
  std::vector<int>         delayed_location_successors;

  std::vector<std::pair<int,int>> local_cyclic_dependencies;

  //======================================== Default constructor
  SPDS()
  {  }

  //======================================== Destructor
//  ~SPDS()
//  {
//    delete spls;
//    for (int gs=0; gs<global_sweep_planes.size(); gs++)
//    {
//      delete global_sweep_planes[gs];
//    }
//  }

  int MapLocJToPrelocI(int locJ)
  {
    for (int i=0; i<location_dependencies.size(); i++)
    {
      if (location_dependencies[i] == locJ)
      {
        return i;
      }
    }

    for (int i=0; i<delayed_location_dependencies.size(); i++)
    {
      if (delayed_location_dependencies[i] == locJ)
      {
        return -(i+1);
      }
    }

    chi_log.Log(LOG_ALLERROR)
      << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
    exit(EXIT_FAILURE);
  }

  int MapLocJToDeplocI(int locJ)
  {
    for (int i=0; i<location_successors.size(); i++)
    {
      if (location_successors[i] == locJ)
      {
        return i;
      }
    }

    chi_log.Log(LOG_ALLERROR)
      << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
    exit(EXIT_FAILURE);
    return -1;
  }


  void AddLocalDependecy(int location_index)
  {
    if (location_index<0){return;}

    bool already_there = false;
    size_t num_deps = location_dependencies.size();
    for (int k=0; k<num_deps; k++)
    {
      if (location_dependencies[k] == location_index)
      {
        already_there = true;
        break;
      }
    }
    if (!already_there)
    {
      location_dependencies.push_back(location_index);
    }
  }

  void AddLocalSuccessor(int location_index)
  {
    if (location_index<0){return;}

    bool already_there = false;
    size_t num_sucs = location_successors.size();
    for (int k=0; k<num_sucs; k++)
    {
      if (location_successors[k] == location_index)
      {
        already_there = true;
        break;
      }
    }
    if (!already_there)
    {
      location_successors.push_back(location_index);
    }
  }
};

#endif