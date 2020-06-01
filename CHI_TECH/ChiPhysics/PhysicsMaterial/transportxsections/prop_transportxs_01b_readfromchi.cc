#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**This method populates a transport cross-section from
 * a Chi cross-section file.*/
void chi_physics::TransportCrossSections::
  MakeFromCHIxsFile(const std::string &file_name)
{

  chi_log.Log(LOG_0) << "Reading Chi cross-section file \"" << file_name << "\"\n";
  //opens and checks if open
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open())
  {
      chi_log.Log(LOG_ALLERROR)<< "Failed to open chi cross-section file \""
          << file_name << "\" in call to "
          << "TransportCrossSections::MakeFromChixsFile\n";
      exit(EXIT_FAILURE);
  }

  //line is used to get rid of lines and word is used to get rid of words
  char line [250];
  std::string word;
  std::string sectionChecker;
  //num groups
  int G = 15;
  //num moments
  int M = 3;

  //ensures the vectors are the correct size
  sigma_tg.resize(G,0.0);
  sigma_fg.resize(G,0.0);
  sigma_captg.resize(G,0.0);
  chi_g.resize(G,0.0);
  nu_sigma_fg.resize(G,0.0);
  
  
  //get rid of the first line of comments
  file.getline(line,250);

  //group and moments
  file >> word >> G;
  file >> word >> M;

  //function for reading in the 1d vectors
  auto Read1DXS = [](std::vector<double>& xs, std::ifstream& file, int G)
  {
      for (int g=0; g<G; g++)
      {
        double group,sigma;
        file>>group>>sigma;
        xs[group] = sigma;
      }
  };

  //reads each section of the 1d xs
  while(sectionChecker!="TRANSFER_MOMENTS_BEGIN")
  {
    file >> sectionChecker;
    if (sectionChecker=="SIGMA_T_BEGIN"){Read1DXS(sigma_tg,file,G);}
    if (sectionChecker=="SIGMA_F_BEGIN"){Read1DXS(sigma_fg,file,G);}
    if (sectionChecker=="NU_BEGIN"){Read1DXS(nu_sigma_fg,file,G);}
  }

  //changes nu_sigma_fg from nu to nu * sigma_fg
  for (int i = 0; i<G;++i){
    nu_sigma_fg[i] = nu_sigma_fg[i]*sigma_fg[i];
  }

  //initilize sparce matrix
  transfer_matrix.resize(M,chi_math::SparseMatrix(G,G));

  while(sectionChecker!="TRANSFER_MOMENTS_END"){
    file>>sectionChecker;
    if(sectionChecker == "TRANSFER_MOMENT_BEGIN"){
      int m,gprime,g;
      file>>m;
      for (int j = 0;j<G;++j){
        file>>word>>gprime>>g;
        for(int i = 0;i<g;++i){
          double destination, value;
          file>>destination>>value;

          transfer_matrix[m].Insert(gprime,g,value);
        }

      }
    }
  }

  file.close();
}