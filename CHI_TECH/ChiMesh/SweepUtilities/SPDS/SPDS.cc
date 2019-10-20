#include "SPDS.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/** Given a location J index, maps to a predecessor location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToPrelocI(int locJ)
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

//###################################################################
/** Given a location J index, maps to a dependent location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToDeplocI(int locJ)
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

//###################################################################
/** Adds a location dependency for this location.*/
void chi_mesh::sweep_management::SPDS::AddLocalDependecy(int location_index)
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

//###################################################################
/** Adds a location successor for this location.*/
void chi_mesh::sweep_management::SPDS::AddLocalSuccessor(int location_index)
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