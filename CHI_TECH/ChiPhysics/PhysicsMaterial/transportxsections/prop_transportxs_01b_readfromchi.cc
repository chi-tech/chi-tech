#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**This method populates a transport cross-section from
 * a PDT cross-section file.*/
void chi_physics::TransportCrossSections::
  MakeFromCHIxsFile(const std::string &file_name)
{

  std::cout << "Reading PDT cross-section file \"" << file_name << "\"\n";
  //opens and checks if open
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open())
  {
      std::cout<< "Failed to open PDT cross-section file \""
          << file_name << "\" in call to "
          << "TransportCrossSections::MakeFromPDTxsFile\n";
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

  
  sigma_tg.resize(G,0.0);
  sigma_fg.resize(G,0.0);
  //sigma_captg.resize(G,0.0); // this is used for SIGMA_S in 1D XS
  chi_g.resize(G,0.0);
  nu_sigma_fg.resize(G,0.0);
  
  
  //get rid of the first line of comments
  file.getline(line,250);
  std::cout<<std::endl<< line <<std::endl;
  //group and moments
  file >> word >> G;
  file >> word >> M;

  auto Read1DXS = [](std::vector<double>& xs, std::ifstream& file, int G)
  {
      for (int g=0; g<G; g++)
      {
        double group,sigma;
        file>>group>>sigma;
        xs[group] = sigma;
      }
  };

  while(sectionChecker!="TRANSFER_MOMENTS_BEGIN")
  {
    file >> sectionChecker;
    if (sectionChecker=="SIGMA_T_BEGIN"){Read1DXS(sigma_tg,file,G);}
    //if (sectionChecker=="SIGMA_S_BEGIN"){Read1DXS(sigma_captg,file,G);}
    if (sectionChecker=="SIGMA_F_BEGIN"){Read1DXS(sigma_fg,file,G);}
    if (sectionChecker=="NU_BEGIN"){Read1DXS(nu_sigma_fg,file,G);}
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
          std::cout<<destination<<" "<<value<<std::endl;
          transfer_matrix[m].Insert(gprime,g,value);
        }
        std::cout<<std::endl;
      }
    }
  }
  

//     //TRANSFER_MOMENT_BEGIN 0
//     int m=0;
     
// // GPRIME_TO_G 2 2
// // 1   0.00157223178365
// // 2   0.365296247644
//     int gprime = 2;
//     //for read 2 lines
//     int g = 1;
//     double val = 0.00157223178365;

//     transfer_matrix[m].Insert(gprime,g,val);

  file.close();
}