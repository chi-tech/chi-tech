#ifndef CHI_SPLS_H
#define CHI_SPLS_H

#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

//###################################################################
/**Contains the sweep plane data.*/
struct chi_mesh::sweep_management::SPLS
{
  std::vector<int> item_id;
};

//###################################################################
/**Stage-wise Task Dependency Graph.
 * Contains the global sweep plane data.*/
struct chi_mesh::sweep_management::STDG
{
  std::vector<int> item_id;
};



#endif //CHI_SPLS_H