#include "fieldfunction.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>

#include <chi_log.h>
#include <chi_mpi.h>
//#include <iostream>
#include <fstream>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

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

/** Writes the VTK "Assembly file" for multiple vtu files.

\param base_filename Base name for all vtu file. .pvtu will get appended to it.
\param field_name Name of the field to be exported
\param num_grps Optional. Defaults to 0. If greater than zero then exports groups.
\author Jason*/
void chi_physics::FieldFunction::WritePVTU(std::string base_filename,
                                           std::string field_name,
                                           int num_grps)
{
    std::string summary_file_name = base_filename + std::string(".pvtu");
    std::ofstream ofile;
    ofile.open(summary_file_name);

    ofile << "<?xml version=\"1.0\"?>" << std::endl;
    ofile << "<!--" << std::endl;
    ofile << "#Unstructured Mesh" << std::endl;
    ofile << "-->" << std::endl;
    ofile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
          << "byte_order=\"LittleEndian\">" << std::endl;
    ofile << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    ofile << "    <PPointData Scalars=\"scalars\">" << std::endl;

    if (num_grps == 0)
    {
        ofile << "      <PDataArray type=\"Float64\" Name=\""
              << field_name << "\" format=\"ascii\"/>" << std::endl;
    }
    else
    {
        for (int g=0; g<num_grps; g++)
        {
            char group_text[100];
            sprintf(group_text,"%03d",g);
            ofile << "      <PDataArray type=\"Float64\" Name=\""
                  << field_name + std::string("_g") + std::string(group_text)
                  << "\" format=\"ascii\"/>" << std::endl;
        }
    }

    ofile << "    </PPointData>" << std::endl;
    ofile << "    <PCellData Scalars=\"scalars\">" << std::endl;
    ofile << "      <PDataArray type=\"Int32\" Name=\"Material\" "
          << " format=\"ascii\"/>" << std::endl;
    ofile << "      <PDataArray type=\"Int32\" Name=\"Partition\""
          << " format=\"ascii\"/>" << std::endl;

    if (num_grps == 0)
    {
        ofile << "      <PDataArray type=\"Float64\" Name=\""
              << field_name + std::string("-Avg") << "\" format=\"ascii\"/>" << std::endl;
    }
    else
    {
        for (int g=0; g<num_grps; g++)
        {
            char group_text[100];
            sprintf(group_text,"%03d",g);
            ofile << "      <PDataArray type=\"Float64\" Name=\""
                  << field_name + std::string("_g") + std::string(group_text) + std::string("_avg")
                  << "\" format=\"ascii\"/>" << std::endl;
        }
    }

    ofile << "    </PCellData>" << std::endl;
    ofile << "    <PPoints>" << std::endl;
    ofile << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << std::endl;
    ofile << "    </PPoints>" << std::endl;

    bool is_global_mesh =
      chi_mesh::GetCurrentHandler()->volume_mesher->options.mesh_global;

    for (int p=0; p<chi_mpi.process_count; p++)
    {
      if (is_global_mesh and p!=0) continue;

      ofile << "      <Piece Source=\""
            << base_filename +
               std::string("_") +
               std::to_string(p) +
               std::string(".vtu")
            << "\"/>" << std::endl;
    }

    ofile << "  </PUnstructuredGrid>" << std::endl;
    ofile << "</VTKFile>" << std::endl;

    ofile.close();
}
