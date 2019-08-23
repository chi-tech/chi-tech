#include "chi_fieldfunction.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern CHI_LOG chi_log;

//###################################################################
/**Exports a field function to VTK format.
 *
 * */
void chi_physics::FieldFunction::ExportToVTK(std::string base_name,
                                             std::string field_name)
{
  chi_log.Log(LOG_0)
    << "Exporting field function " << text_name
    << " to files with base name " << base_name;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PWLD NODES
  if (type == FF_SDM_FV)
    ExportToVTKFV(base_name,field_name);
  if (type == FF_SDM_CFEM)
    ExportToVTKPWLC(base_name,field_name);
  if (type == FF_SDM_PWLD)
    ExportToVTKPWLD(base_name,field_name);

}

//###################################################################
/**Exports a field function to VTK format but exports all of the
 * available groups.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKG(std::string base_name,
                                             std::string field_name)
{
  chi_log.Log(LOG_0)
    << "Exporting field function " << text_name
    << " to files with base name " << base_name;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PWLD NODES
  if (type == FF_SDM_FV)
    ExportToVTKFVG(base_name,field_name);
  if (type == FF_SDM_CFEM)
    ExportToVTKPWLCG(base_name,field_name);
  if (type == FF_SDM_PWLD)
    ExportToVTKPWLDG(base_name,field_name);

}
