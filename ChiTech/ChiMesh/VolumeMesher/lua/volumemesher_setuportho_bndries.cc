#include "ChiLua/chi_lua.h"

#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

//###################################################################
/** Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as xmax=0, xmin=1, ymax=2, ymin=3, zmax=4, zmin=5.*/
int chiVolumeMesherSetupOrthogonalBoundaries(lua_State* L)
{
  chi_mesh::VolumeMesher::SetupOrthogonalBoundaries();
  return 0;
}