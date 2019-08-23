#include "CHI_PHYSICS/CHI_PHYSICSMATERIAL/property10_transportxsections.h"

//#include <unistd.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern CHI_LOG chi_log;
extern CHI_MPI chi_mpi;

//###################################################################
/**This method populates a transport cross-section from
 * a PDT cross-section file.*/
void chi_physics::TransportCrossSections::
  MakeFromPDTxsFile(const std::string &file_name,std::string MT_TRANSFER)
{
  chi_log.Log(LOG_0)
    << "Reading PDT cross-section file \"" << file_name << "\"";

  //===================================================== Opening the file
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
  transfer_matrix.clear();
  transfer_matrix.resize(scat_order+1,
                         chi_math::SparseMatrix(num_grps_G,num_grps_G));

  //====================================== Skip ahead to fine "MT 1"
  std::string word0,word1;
  std::stringstream linestring;

  while (!file.eof())
  {
    file.getline(line,250);
    linestring = std::stringstream(line);

    linestring >> word0 >> word1;

    if ((word0 == "MT") && (word1 == "1")) break;

    if (file.eof())
    {
      chi_log.Log(LOG_0ERROR)
        << "While reading file \"" << file_name << "\""
        << " premature end-of-file while searching for MT 1"
        << " or MT 1 was not found";
      exit(EXIT_FAILURE);
    }
  }

  for (int g=0; g<G; g++)
  {
    double value=0.0;
    file >> value;

    sigma_tg[g] = value;
  }


  //======================================== Find all moments of transfer matrix
  std::string MT_SEARCH_VAL = MT_TRANSFER + std::string(",");
  int mom=-1;
  while (!file.eof())
  {
    file.getline(line,250);
    linestring = std::stringstream(line);

    linestring >> word0 >> word1;


    if ((word0 == "MT") && (word1 == MT_SEARCH_VAL))
    {
      word0="";
      word1="";
      mom++;
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

//  chi_log.Log(LOG_0)
//  << "Transfer matrix mom 0 group 167:\n"
//  << transfer_matrix[0].PrintS();



  file.close();
}