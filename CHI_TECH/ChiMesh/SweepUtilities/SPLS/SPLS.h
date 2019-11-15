#ifndef _chi_spls_h
#define _chi_spls_h

#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

//###################################################################
/**Contains the sweep plane data.*/
struct chi_mesh::sweep_management::SPLS
{
  std::vector<int> item_id;
  //chi_mesh::sweep_management::FLUDS* fluds;
};

//###################################################################
/**Stage-wise Task Dependency Graph.
 * Contains the global sweep plane data.*/
struct chi_mesh::sweep_management::STDG
{
  std::vector<int> item_id;

};



#endif