#include "surfacemesher.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Virtual execute function. Meant to be overwritten.*/
void chi_mesh::SurfaceMesher::Execute()
{
  chi::log.Log0Verbose1() << "This is an empty mesher. Nothing to execute.";
}