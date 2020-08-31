#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"

#include <chi_log.h>

extern ChiLog& chi_log;

#include <string>

/**\defgroup ChiXSFile Chi-Tech Cross-section format
 *\ingroup LuaPhysicsMaterials
 *
 * An example Chi-Tech cross-section file is shown below:
\code
# Bla bla bla
# This header can be as large as you please. The actual processing
# starts at NUM_GROUPS as the first word. After that, NUM_MOMENTS needs to
# be processed before any of the other keywords.
NUM_GROUPS 2
NUM_MOMENTS 2
SIGMA_T_BEGIN
0   0.5
1   0.5
SIGMA_T_END

Comments Bla Bla


SIGMA_F_BEGIN
0   0.01
1   0.40737
SIGMA_F_END

NU_BEGIN
0    2.5
1    2.5
NU_END
CHI_BEGIN
0    1.0
1    0.0
CHI_END

# 1/v terms
DDT_COEFF_BEGIN
0    4.5454e-10   #1/2.2e10
1    3.6745e-3
DDT_COEFF_END


TRANSFER_MOMENTS_BEGIN
#Zeroth moment (l=0)
M_GPRIME_G_VAL 0 0 0 0.01
M_GPRIME_G_VAL 0 0 1 0.09
M_GPRIME_G_VAL 0 1 1 0.08

#(l=1)
M_GPRIME_G_VAL 1 0 0 -0.001
M_GPRIME_G_VAL 1 0 1 0.001
M_GPRIME_G_VAL 1 1 1 0.001
TRANSFER_MOMENTS_END
\endcode

## Details:

- NUM_GROUPS     G       Required. Specified the number of groups for this cross-section
- NUM_MOMENTS   M    Required. Transfer matrices to allocate for this cross-section (whether used or not).
- Optional key words per line:
  - SIGMA_T_BEGIN.   Starts a block that is terminated by a line SIGMA_T_END. Each line in the block processes the first two words as [Group, sigma_t]. Populates the sigma_tg field.
  - SIGMA_F_BEGIN.   Starts a block that is terminated by a line SIGMA_F_END. Each line in the block processes the first two words as [Group, sigma_f]. Populates the sigma_fg field.
  - NU_BEGIN.   Starts a block that is terminated by a line NU_END. Each line in the block processes the first two words as [Group, nu]. Populates the nu_sigma_fg field. Upon completing the file processing the nu_sigma_fg field gets multiplied by sigma_fg.
  - CHI_BEGIN.   Starts a block that is terminated by a line CHI_END. Each line in the block processes the first two words as [Group, chi]. Populates the chi_g field.
  - DDT_COEFF_BEGIN.   Starts a block that is terminated by a line DDT_COEFF_END. Each line in the block processes the first two words as [Group, ddt_coeff]. Populates the ddt_coeff field.
  - TRANSFER_MOMENTS_BEGIN. Starts a block that is terminated by a line TRANSFER_MOMENTS_END. Each line in the block processes a line only if it starts with the keyword M_GPRIME_G_VAL which needs to be followed by four values [moment,gprime,g,value]. Populates transfer-matrix for moment m, row g, column gprime, with value.
- Comments can be between individual blocks but only the TRANSFER_MOMENTS block may have comments between the _BEGIN and _END
- Comments do not have to start with any specific character since the file format is keyword driven.
- Any number that is not convertible to its required type (integer, double) will throw an error to that effect.
- All errors will indicate the current file, line number and nature of the error.

 * */

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
  std::string line;
  std::string word, first_word;
  std::string sectionChecker;

  //num moments
  int M = L;

  //#############################################
  /**Lambda for converting strings to doubles.*/
  auto StrToD = [](const char* str)
  {
    char* endptr;
    double value = std::strtod(str, &endptr);

    if (endptr[0] != 0)
      throw std::runtime_error(
        std::string("Invalid floating-point "
                    "number conversion: \"")+str+"\".");

    return value;
  };

  //#############################################
  /**Lambda for converting strings to integers.*/
  auto StrToI = [](const char* str)
  {
    char* endptr;
    long int value = std::strtol(str, &endptr,10);

    if (endptr[0] != 0)
      throw std::runtime_error(
        std::string("Invalid integer "
                    "number conversion: \"")+str+"\".");

    return value;
  };

  //#############################################
  /**Lambda function for reading in the 1d vectors.*/
  auto Read1DXS = [StrToD,StrToI]
    (std::string keyword,std::vector<double>& xs,
     std::ifstream& file, int Gtot, int& line_number,
     std::istringstream& line_stream)
  {
    int g=-1;
    char first_word[250];
    char line[250];
    char value_str[250];

    file.getline(line,250); ++g; ++line_number;
    line_stream = std::istringstream(line);
    line_stream >> first_word;

    while (std::string(first_word) != (keyword + "_END"))
    {
      if (g>=Gtot) throw std::runtime_error(
        "Too many lines in block " + keyword + ".");

      line_stream >> value_str;

      int group    = StrToI(first_word);
      double value = StrToD(value_str);

      xs[group] = value;

      file.getline(line,250); ++g; ++line_number;
      line_stream = std::istringstream(line);
      line_stream >> first_word;
    }
  };

  //#############################################
  /**Lambda reading a transfer matrix.*/
  auto ReadTransferMatrix = [StrToD,StrToI]
    (std::string keyword,std::vector<chi_math::SparseMatrix>& matrix,
     std::ifstream& file, int Gtot, int& line_number,
     std::istringstream& line_stream)
  {
    char first_word[250];
    char line[250];
    char value_str0[250],value_str1[250],value_str2[250],value_str3[250];

    file.getline(line,250); ++line_number;
    line_stream = std::istringstream(line);
    line_stream >> first_word;

    while (std::string(first_word) != (keyword + "_END") and !file.eof())
    {
      if (std::string(first_word) == "M_GPRIME_G_VAL")
      {
        value_str0[0] = '\0';
        value_str1[0] = '\0';
        value_str2[0] = '\0';
        value_str3[0] = '\0';
        line_stream >> value_str0 >> value_str1 >> value_str2 >> value_str3;

        if (value_str0[0]==0 or value_str1[0]==0 or
            value_str2[0]==0 or value_str3[0]==0)
          throw std::runtime_error("Invalid amount of arguments. "
                                   "Requires 4 numbers.");

        int m        = StrToI(value_str0);
        int gprime   = StrToI(value_str1);
        int g        = StrToI(value_str2);
        double value = StrToD(value_str3);

        if (m<matrix.size())
          matrix[m].Insert(g,gprime,value);
      }

      file.getline(line,250); ++line_number;
      line_stream = std::istringstream(line);
      line_stream >> first_word;
    }
  };

  //================================================== Read file line by line
  bool grabbed_G = false;
  int line_number=0;
  bool not_eof = bool(std::getline(file,line)); ++line_number;
  while (not_eof)
  {
    std::istringstream line_stream(line);
    line_stream >> first_word;

//    chi_log.Log(LOG_0) << line << " _" << first_word << "_";

    if (first_word == "NUM_GROUPS") {line_stream >> G; grabbed_G = true;}
    if (first_word == "NUM_MOMENTS")
    {
      line_stream >> M;
      if (grabbed_G)
      {
        sigma_tg.resize(G,0.0);
        sigma_fg = sigma_captg = chi_g = nu_sigma_fg = ddt_coeff = sigma_tg;
        transfer_matrix.resize(M,chi_math::SparseMatrix(G,G));
      }
    }

    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f = file;
      auto& fw = first_word;

      if (fw == "SIGMA_T_BEGIN")   Read1DXS ("SIGMA_T",sigma_tg,f,G,ln,ls);
      if (fw == "SIGMA_F_BEGIN")   Read1DXS("SIGMA_F",sigma_fg,f,G,ln,ls);
      if (fw == "NU_BEGIN")        Read1DXS("NU",nu_sigma_fg,f,G,ln,ls);
      if (fw == "CHI_BEGIN")       Read1DXS("CHI",chi_g,f,G,ln,ls);
      if (fw == "DDT_COEFF_BEGIN") Read1DXS("DDT_COEFF",ddt_coeff,f,G,ln,ls);

      if (fw == "TRANSFER_MOMENTS_BEGIN")
        ReadTransferMatrix("TRANSFER_MOMENTS",transfer_matrix,f,G,ln,ls);

    }

    catch (const std::runtime_error& err)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Error reading xs-file \"" + file_name + "\". "
        << "Line number " << line_number << ". "
        << err.what();
      exit(EXIT_FAILURE);
    }
    catch (...)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Nope";
      exit(EXIT_FAILURE);
    }

    first_word = "";
    not_eof = bool(std::getline(file,line)); ++line_number;
  }
  L = M-1;

  //changes nu_sigma_fg from nu to nu * sigma_fg
  for (int i = 0; i<G;++i){
    nu_sigma_fg[i] = nu_sigma_fg[i]*sigma_fg[i];
  }

  file.close();
}