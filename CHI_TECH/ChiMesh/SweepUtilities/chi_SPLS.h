#ifndef _chi_spls_h
#define _chi_spls_h

#include "chi_FLUDS.h"

//###################################################################
/**Contains the sweep plane data.*/
struct chi_mesh::SweepManagement::SPLS
{
  std::vector<int> item_id;
  //chi_mesh::SweepManagement::FLUDS* fluds;
};

//###################################################################
/**Stage-wise Task Dependency Graph.
 * Contains the global sweep plane data.*/
struct chi_mesh::SweepManagement::STDG
{
  std::vector<int> item_id;

};



#endif