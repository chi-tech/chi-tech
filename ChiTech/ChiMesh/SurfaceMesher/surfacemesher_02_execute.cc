#include "surfacemesher.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Virtual execute function. Meant to be overwritten.*/
void chi_mesh::SurfaceMesher::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "This is an empty mesher. Nothing to execute.";
}