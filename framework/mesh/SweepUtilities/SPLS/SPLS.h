#ifndef CHI_SPLS_H
#define CHI_SPLS_H

#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

namespace chi_mesh::sweep_management
{

// ###################################################################
/**Contains the sweep plane data.*/
struct SPLS
{
  std::vector<int> item_id;
};

// ###################################################################
/**Stage-wise Task Dependency Graph.
 * Contains the global sweep plane data.*/
struct STDG
{
  std::vector<int> item_id;
};

} // namespace chi_mesh::sweep_management

#endif // CHI_SPLS_H