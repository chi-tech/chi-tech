#include "material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <string>


/**\defgroup ChiXSFile Chi-Tech Cross-section format 1
 *\ingroup LuaPhysicsMaterials
 *
 * An example Chi-Tech cross-section file is shown below. The bare-bones
 * format is shown below with more examples below:
\code
# This header can be as large as you please. The actual processing
# starts at NUM_GROUPS as the first word. After that, NUM_MOMENTS needs to
# be processed before any of the other keywords.
NUM_GROUPS 2
NUM_MOMENTS 2
SIGMA_T_BEGIN
0   0.5
1   0.5
SIGMA_T_END

Comments

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

## Steady state simulations:

The cross sections can be used in simulations solving the steady state
Linear Boltzmann Equation of the form:

\f[
\vec{\Omega}_n \boldsymbol{\cdot} \vec{\nabla} \psi_{ng} +
\sigma_{tg} \psi_{ng} =
\sum_{\ell=0}^L \sum_{m=-\ell}^{+\ell} \frac{2\ell+1}{4\pi}
    Y_{\ell m} (\vec{\Omega}_n)
    \sum_{g'=0}^{G-1} \sigma_{s\ell,g'{\to}g} \phi_{\ell m,g'} +
    q_{ext,ng} + q_{fission,ng}
\f]

The two most prominent items required here includes \f$ \sigma_{tg} \f$ and
 \f$ \sigma_{s\ell,g'{\to}g} \f$. The latter is an entry in a structure we
 call a generic transfer matrix for moment \f$ \ell \f$ with rows \f$g=0..G-1\f$
 and columns \f$g'=0..G-1\f$. These two items are often the only items required
 in a transport simulation.

 In simulations with fission-sources, cross sections
 support two formats, the simple combined cross sections without delayed
 neutrons

\f[
q_{fission,ng} = \frac{\chi_g}{4\pi}
 \sum_{g'=0}^{G-1} \nu_{g'} \sigma_{fg'} \phi_{00g'}
\f]

and those with delayed neutrons, the latter which are currently only used in the
KEigen solver.

## KEigen-value related items
As stated before, the cross section file supports two formats, the simple combined cross sections without delayed
neutrons which are depicted above and those with delayed neutrons shown below

\f[
q_{fission,ng} = \frac{\chi_g}{4\pi}
 \sum_{g'=0}^{G-1} \nu_{prompt,g'} \sigma_{fg'} \phi_{00g'}
 + \sum_{j=0}^{J-1} \frac{\chi_{delayed,jg}}{4\pi} \gamma_j
    \sum_{g'=0}^{G-1} \nu_{delayed,g'} \sigma_{fg'} \phi_{00g'}
\f].

Codes may also choose to update precursor concentrations for which their decay
constants, \f$ \lambda_j \f$ are required.

## Keyword definitions

- NUM_GROUPS num_groups Required. Specifies the number of groups for this
  cross-section. Symbol \f$ G \f$.
- NUM_MOMENTS num_moments Required. The number of transfer matrices to allocate
  for this cross-section (whether used or not). Typically this number is one
  greater than the scattering order (i.e., \f$ L+1 \f$)
- NUM_PRECURSORS num_precursors Optional. Indicates how many precursors are used
  in this cross section. Symbol \f$ J \f$
- Optional key words per line:
  - SIGMA_T_BEGIN. Optional. Starts a block that is terminated by a line SIGMA_T_END.
    Each line in the block processes the first two words as [Group, sigma_t].
    Populates the sigma_tg field. Symbol \f$ \sigma_{tg} \f$.
  - SIGMA_F_BEGIN. Optional. Starts a block that is terminated by a line SIGMA_F_END.
    Each line in the block processes the first two words as [Group, sigma_f].
    Populates the sigma_fg field. Symbol \f$ \sigma_{fg} \f$.
  - SIGMA_A_BEGIN. Optional. Starts a block that is terminated by a line SIGMA_A_END.
    Each line in the block processes the first two words as [Group, sigma_a].
    Populates the sigma_ag field. Symbol \f$ \sigma_{ag} \f$. If this is
    not supplied then sigma_a is estimated from the transfer matrix and
    may erroneously estimate balance.
  - NU_BEGIN. Optional. Starts a block that is terminated by a line NU_END. Each line in
    the block processes the first two words as [Group, nu]. Populates the
    nu field. Upon completing the file processing the nu_sigma_fg field
    gets populated from the product of nu and sigma_fg. Symbol \f$ \nu_g \f$.
  - NU_PROMPT_BEGIN. Optional. Starts a block that is terminated by a line
    NU_PROMPT_END. Each line in the block processes the first two words as
    [Group, nu]. Populates the nu_prompt field. Upon completing the file
    processing the nu_p_sigma_fg field gets populated from the product of
    nu_prompt and sigma_fg.
    Symbol \f$ \nu_{prompt,g} \f$.
  - NU_DELAYED_BEGIN. Optional. Starts a block that is terminated by a line
    NU_DELAYED_END. Each line in the block processes the first two words as
    [Group, nu]. Populates the nu_delayed field. Upon completing the file
    processing the nu_d_sigma_fg field gets populated from the product of
    nu_delayed and sigma_fg.
    Symbol \f$ \nu_{delayed,g} \f$.
  - CHI_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_END. Each line in the block processes the first two words as
    [Group, chi]. Populates the chi field. Symbol \f$ \chi_{g} \f$.
  - CHI_PROMPT_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_PROMPT_END. Each line in the block processes the first two words as
    [Group, chi]. Populates the chi_prompt field. Symbol \f$ \chi_{prompt, g} \f$.
  - VELOCITY_BEGIN. Optional. Starts a block that is terminated by a line
    VELOCITY_END. Each line in the block processes the first two words as
    [Group, velocity]. Populates the velocity field.
    Symbol \f$ \frac{1}{v_g} \f$.

  - PRECURSOR_LAMBDA_BEGIN. Optional. Starts a block that is terminated by a line
    PRECURSOR_LAMBDA_END. Each line in the block processes the first two words as
    [precursor, lambda]. Populates the lambda field (the precursor decay
    constant). Symbol \f$ \lambda_j \f$.
  - PRECURSOR_YIELD_BEGIN. Optional. Starts a block that is terminated by a line
    PRECURSOR_YIELD_END. Each line in the block processes the first two words as
    [precursor, gamma]. Populates the gamma field (the precursor production
    fraction per fission). Symbol \f$ \gamma_j \f$.

  - CHI_DELAYED_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_DELAYED_END. Each line in the block processes the first word as the
    group index and the remaining NUM_PRECURSORS words as the the individual
    precursor's associated delayed spectrum (chi). Populates the chi_d field.
    Symbol \f$ \chi_{delayed,jg} \f$.

  - TRANSFER_MOMENTS_BEGIN. Optional. Starts a block that is terminated by a line
    TRANSFER_MOMENTS_END. Each line in the block processes a line only if it
    starts with the keyword M_GPRIME_G_VAL which needs to be followed by four
    values [moment,gprime,g,value]. Populates transfer-matrix for moment m,
    row g, column gprime, with value. Symbol \f$ \sigma_{s\ell,g'{\to}g} \f$.
- Comments can be between individual blocks but only the TRANSFER_MOMENTS block
  may have comments between the _BEGIN and _END
- Comments do not have to start with any specific character since the file
  format is keyword driven.
- Any number that is not convertible to its required type (integer, double)
  will throw an error to that effect.
- All errors will indicate the current file, line number and nature of the error.

## More Advanced Examples
\code
# This header can be as large as you please. The actual processing
# starts at NUM_GROUPS as the first word. After that, NUM_MOMENTS needs to
# be processed before any of the other keywords.
NUM_GROUPS 2
NUM_MOMENTS 2
NUM_PRECURSORS 3
SIGMA_T_BEGIN
0   0.5
1   0.5
SIGMA_T_END

Comments

SIGMA_F_BEGIN
0   0.01
1   0.40737
SIGMA_F_END

NU_PROMPT_BEGIN
0    2.45
1    2.45
NU_PROMPT_END

CHI_PROMPT_BEGIN
0    1.0
1    0.0
CHI_PROMPT_END

VELOCITY_BEGIN
0    2.2e10
1    272.145
VELOCITY_END

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

PRECURSOR_LAMBDA_BEGIN
0		0.1
1   0.2
2   0.3
PRECURSOR_LAMBDA_END

PRECURSOR_GAMMA_BEGIN
0		0.25
1   0.5
2   0.25
PRECURSOR_GAMMA_END

NU_DELAYED_BEGIN
0		0.01
1   0.02
2   0.01
NU_DELAYED_END

CHI_DELAYED_BEGIN
G_PRECURSORJ_VAL 0  0	1.0
G_PRECURSORJ_VAL 0  1	1.0
G_PRECURSORJ_VAL 0  2	1.0

G_PRECURSORJ_VAL 1  0	0.0
G_PRECURSORJ_VAL 1  1	0.0
G_PRECURSORJ_VAL 1  2	0.0
CHI_DELAYED_END
\endcode
 * */

//###################################################################
/**This method populates a transport cross-section from
 * a Chi cross-section file.*/
void chi_physics::TransportCrossSections::
  MakeFromCHIxsFile(const std::string &file_name)
{
  //======================================== Clear any previous data
  Reset();

  //======================================== Read file
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
  bool grabbed_G = false;
  bool grabbed_M = false;

  //num moments
  size_t M = scattering_order+1; //just to init M

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

  auto ThrowGandMError = [file_name]()
  {
    throw std::runtime_error("ChiTech format cross-section file "
                             + file_name + "is trying to write data without "
                             "first processing NUM_GROUPS and NUM_MOMENTS.");
  };

  //#############################################
  /**Lambda function for reading in the 1d vectors.*/
  auto Read1DXS = [StrToD,StrToI,&grabbed_G,&grabbed_M,ThrowGandMError]
    (const std::string& keyword,std::vector<double>& xs,
     std::ifstream& file, size_t Gtot, int& line_number,
     std::istringstream& line_stream)
  {
    if ((not grabbed_G) or (not grabbed_M)) ThrowGandMError();
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
  /**Lambda reading the chi delayed matrix.*/
  auto ReadChiDelayed = [StrToD,StrToI,&grabbed_G,&grabbed_M,ThrowGandMError]
    (const std::string& keyword,std::vector<std::vector<double>>& spectra,
     std::ifstream& file, size_t Gtot, int& line_number,
     std::istringstream& line_stream)
  {
    if ((not grabbed_G) or (not grabbed_M)) ThrowGandMError();
    char first_word[250];
    char line[250];
    char value_str0[250],value_str1[250],value_str2[250];

    file.getline(line,250); ++line_number;
    line_stream = std::istringstream(line);
    line_stream >> first_word;

    while (std::string(first_word) != (keyword + "_END"))
    {
      if (std::string(first_word) == "G_PRECURSORJ_VAL")
      {
        value_str0[0] = '\0';
        value_str1[0] = '\0';
        value_str2[0] = '\0';
        line_stream >> value_str0 >> value_str1 >> value_str2;

        if (value_str0[0]==0 or value_str1[0]==0 or value_str2[0]==0)
          throw std::runtime_error("Invalid amount of arguments. "
                                   "Requires 4 numbers.");

        int group     = StrToI(value_str0);
        int precursor = StrToI(value_str1);
        double value  = StrToD(value_str2);

        spectra[group][precursor] = value;
      }

      file.getline(line,250); ++line_number;
      line_stream = std::istringstream(line);
      line_stream >> first_word;
      if (line_stream.str().empty())
        sprintf(first_word," ");
    }
  };

  //#############################################
  /**Lambda reading a transfer matrix.*/
  auto ReadTransferMatrix = [StrToD,StrToI,&grabbed_G,&grabbed_M,ThrowGandMError]
    (const std::string& keyword,std::vector<chi_math::SparseMatrix>& matrix,
     std::ifstream& file, size_t Gtot, int& line_number,
     std::istringstream& line_stream)
  {
    if ((not grabbed_G) or (not grabbed_M)) ThrowGandMError();
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
                                   "Requires 4 numbers.\nLine:" + line_stream.str());

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
      if (line_stream.str().empty())
        sprintf(first_word," ");
    }
  };

  //================================================== Read file line by line
  int line_number=0;
  bool not_eof = bool(std::getline(file,line)); ++line_number;
  while (not_eof)
  {
    std::istringstream line_stream(line);
    line_stream >> first_word;

    if (first_word == "NUM_GROUPS")
    {
      line_stream >> num_groups;
      grabbed_G = true;
      sigma_t.assign(num_groups, 0.0);
      sigma_f.assign(num_groups, 0.0);
      sigma_a.assign(num_groups, 0.0);
      chi.assign(num_groups, 0.0);
      chi_prompt.assign(num_groups, 0.0);
      nu.assign(num_groups, 0.0);
      nu_prompt.assign(num_groups, 0.0);
      nu_delayed.assign(num_groups, 0.0);
      nu_sigma_f.assign(num_groups, 0.0);
      nu_prompt_sigma_f.assign(num_groups, 0.0);
      nu_delayed_sigma_f.assign(num_groups, 0.0);
      velocity.assign(num_groups, 0.0);
    }
    if (first_word == "NUM_MOMENTS")
    {
      line_stream >> M;
      grabbed_M = true;
      if (grabbed_G)
        transfer_matrices.resize(M, chi_math::SparseMatrix(num_groups, num_groups));
    }
    if (first_word == "NUM_PRECURSORS")
    {
      line_stream >> num_precursors;
      precursor_lambda.resize(num_precursors, 0.0);
      precursor_yield.resize(num_precursors, 0.0);

      if (grabbed_G)
      {
        chi_delayed.resize(num_groups);
        for (int g=0; g < num_groups; ++g)
          chi_delayed[g].resize(num_precursors);
      }
    }

    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f  = file;
      auto& fw = first_word;

      if (fw == "SIGMA_T_BEGIN")       Read1DXS ("SIGMA_T"    , sigma_t    , f, num_groups, ln, ls);
      if (fw == "SIGMA_F_BEGIN")       Read1DXS("SIGMA_F"     , sigma_f    , f, num_groups, ln, ls);
      if (fw == "SIGMA_A_BEGIN")       Read1DXS("SIGMA_A"     , sigma_a    , f, num_groups, ln, ls);
      if (fw == "NU_BEGIN")            Read1DXS("NU"          , nu         , f, num_groups, ln, ls);
      if (fw == "NU_PROMPT_BEGIN")     Read1DXS("NU_PROMPT"   , nu_prompt  , f, num_groups, ln, ls);
      if (fw == "NU_DELAYED_BEGIN")    Read1DXS("NU_DELAYED"  , nu_delayed , f, num_groups, ln, ls);
      if (fw == "CHI_BEGIN")           Read1DXS("CHI"         , chi        , f, num_groups, ln, ls);
      if (fw == "CHI_PROMPT_BEGIN")    Read1DXS("CHI_PROMPT"  , chi_prompt , f, num_groups, ln, ls);
      if (fw == "VELOCITY_BEGIN")      Read1DXS("VELOCITY"    , velocity   , f, num_groups, ln, ls);

      if (fw == "TRANSFER_MOMENTS_BEGIN")
        ReadTransferMatrix("TRANSFER_MOMENTS",
                           transfer_matrices, f, num_groups, ln, ls);

      if (num_precursors > 0)
      {
        if (fw == "PRECURSOR_LAMBDA_BEGIN")
          Read1DXS("PRECURSOR_LAMBDA", precursor_lambda, f, num_precursors, ln, ls);
        if (fw == "PRECURSOR_YIELD_BEGIN")
          Read1DXS("PRECURSOR_YIELD", precursor_yield, f, num_precursors, ln, ls);
        if (fw == "CHI_DELAYED_BEGIN")
          ReadChiDelayed("CHI_DELAYED", chi_delayed, f, num_groups, ln, ls);
      }
    }//try
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
        << "Unknown error in " << std::string(__FUNCTION__);
      exit(EXIT_FAILURE);
    }

    first_word = "";
    not_eof = bool(std::getline(file,line)); ++line_number;
  }//while not EOF, read each lines
  scattering_order = M-1;

  //perform checks and enforce physical relationships
  ValidateCrossSections();

  //compute nu_sigma_f terms
  for (size_t g = 0; g < num_groups; ++g)
  {
    nu_sigma_f        [g] = nu        [g] * sigma_f[g];
    nu_prompt_sigma_f [g] = nu_prompt [g] * sigma_f[g];
    nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
  }

  file.close();
}