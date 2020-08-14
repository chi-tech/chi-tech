#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**This method populates a transport cross-section from
 * a PDT cross-section file.*/
void chi_physics::TransportCrossSections::
  MakeFromPDTxsFile(const std::string &file_name,std::string MT_TRANSFER)
{
  chi_log.Log(LOG_0)
    << "Reading PDT cross-section file \"" << file_name << "\"";

  std::string MT_SEARCH_VAL = MT_TRANSFER + std::string(",");

  //======================================== Opening the file
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open())
  {
    chi_log.Log(LOG_0ERROR)
      << "Failed to open PDT cross-section file \""
      << file_name << "\" in call to "
      << "TransportCrossSections::MakeFromPDTxsFile";
    exit(EXIT_FAILURE);
  }

  char line[250];
  std::string word;
  std::string mg_or_single_group;
  std::string xs_type;
  int         num_grps_G = 0;
  int         num_process= 0;
  int         num_transfer=0;
  int         scat_order  =0;

  //======================================== Read header information
  //The first two lines are comments
  file.getline(line,250);
  file.getline(line,250);

  //This file is a "multigroup" "neutron" library generated from Python.
  file >> word >> word >> word >> word;
  file >> mg_or_single_group;
  file >> xs_type;
  file.getline(line,250);

  if ((mg_or_single_group != "multigroup") && (xs_type != "neutron"))
  {
    chi_log.Log(LOG_0ERROR)
      << "Currently only multigroup neutron cross-section can be"
      << " read from PDT cross-section files. The file \""
      << file_name << "\" has " << mg_or_single_group << " "
      << xs_type << " cross-sections.";
    file.close();
    exit(EXIT_FAILURE);
  }

  //1 temperatures, 1 densities, and "168" groups.
  file >> word >> word >> word >> word >> word >> num_grps_G;
  file.getline(line,250);
  file.getline(line,250);


  //"6" neutron processes and "1" transfer process.
  file >> num_process >> word >> word >> word >> num_transfer;
  file.getline(line,250);

  //Scattering order "3"
  file >> word >> word >> scat_order;
  file.getline(line,250);

  //====================================== Resizing cross-sections
  G = num_grps_G;
  L = scat_order;
  sigma_tg.clear();
  sigma_tg.resize(num_grps_G,0.0);
  sigma_fg.resize(num_grps_G,0.0);
  sigma_captg.resize(num_grps_G,0.0);
  chi_g.resize(num_grps_G,0.0);
  nu_sigma_fg.resize(num_grps_G,0.0);
  ddt_coeff.resize(num_grps_G, 0.0);

  transfer_matrix.clear();
  transfer_matrix.resize(scat_order+1,
                         chi_math::SparseMatrix(num_grps_G,num_grps_G));

  //======================================== Lambda for advancing to MT
  auto AdvanceToNextMT = [](std::ifstream& file,
                            const std::string& file_name)
  {
    while (!file.eof())
    {
      char line[250];
      file.getline(line,250);
      auto linestring = std::stringstream(line);
      std::string word0,word1;

      linestring >> word0 >> word1;

      if (word0 == "MT")
      {
        int word1_val = std::stoi(word1,nullptr);
        return word1_val;
      }
    }
    return 0;
  };

  //======================================== Lambda for reading 1D xs
  auto Read1DXS = [](std::vector<double>& xs, std::ifstream& file, int G)
  {
    for (int g=0; g<G; g++)
    {
      double value=0.0;
      file >> value;

      xs[g] = value;
    }
  };

  //======================================== Find MT-by-MT
  std::string word0,word1;
  std::stringstream linestring;
  int mt_transfer = std::stoi(MT_SEARCH_VAL, nullptr);
  int mt_number = 1;
  int mom = -1;

  while (mt_number != 0)
  {
    mt_number = AdvanceToNextMT(file,file_name);

    if (mt_number == 1)    Read1DXS(sigma_tg,file,G);
    if (mt_number == 18)   Read1DXS(sigma_fg,file,G);
    if (mt_number == 27)   Read1DXS(sigma_captg,file,G);
    if (mt_number == 2018) Read1DXS(chi_g,file,G);
    if (mt_number == 2452) Read1DXS(nu_sigma_fg,file,G);

    if (mt_number == mt_transfer)
    {
      ++mom;
      for (int g=0; g<G; g++)
      {
        int sink=-1,gprime_first=-1, gprime_last=-2;
        file.getline(line,250);
        linestring = std::stringstream(line);

        //      Sink,   first,  last:   165     104             167
        linestring >> word >> word >> word >> sink >> gprime_first >> gprime_last;

        if (sink != g)
        {
          chi_log.Log(LOG_0ERROR)
            << "Mismatched sink group with general group structure "
               "encountered during transfer moment processing. " << sink
            << " " << g;
          exit(EXIT_FAILURE);
        }

        for (int gprime=gprime_first; gprime<=gprime_last; gprime++)
        {
          double value = 0.0;
          file >> value;
          transfer_matrix[mom].Insert(g,gprime,value);
        }
        file.getline(line,250);
      }
    }
  }


  file.close();
}