#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <string>
#include <algorithm>


/**\defgroup ChiXSFile Chi-Tech Cross-section format 1
 *\ingroup LuaPhysics
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

In simulations with fission-sources, cross sections support two formats,
the simple combined cross sections without delayed neutrons

\f[
q_{fission,ng} = \frac{\chi_g}{4\pi}
 \sum_{g'=0}^{G-1} \nu_{g'} \sigma_{fg'} \phi_{00g'}
\f]

and those with delayed neutrons, the latter which are currently only used in
the k-eigenvalue solver.

## k-eigenvalue related items
As stated before, the cross section file supports two formats, the simple
combined cross sections without delayed neutrons which are depicted above
and those with delayed neutrons shown below

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
  - SIGMA_T_BEGIN. Optional. Starts a block that is terminated by a line
    SIGMA_T_END. Each line in the block processes the first two words as
    [group, sigma_t]. Populates the sigma_tg field. Symbol \f$ \sigma_{tg} \f$.
  - SIGMA_A_BEGIN. Optional. Starts a block that is terminated by a line
    SIGMA_A_END. Each line in the block processes the first two words as
    [group, sigma_a]. Populates the sigma_ag field. Symbol \f$ \sigma_{ag} \f$.
    If this is not supplied then sigma_a is estimated from the transfer matrix
    and may erroneously estimate balance.
  - SIGMA_F_BEGIN. Optional. Starts a block that is terminated by a line
    SIGMA_F_END. Each line in the block processes the first two words as
    [group, sigma_f]. Populates the sigma_fg field. Symbol \f$ \sigma_{fg} \f$.
  - NU_BEGIN. Optional. Starts a block that is terminated by a line NU_END.
    Each line in the block processes the first two words as [group, nu].
    Populates the nu field. Upon completing the file processing the nu_sigma_fg
    field gets populated from the product of nu and sigma_fg.
    Symbol \f$ \nu_g \f$.
  - NU_PROMPT_BEGIN. Optional. Starts a block that is terminated by a line
    NU_PROMPT_END. Each line in the block processes the first two words as
    [group, nu_prompt]. Populates the nu_prompt field. Upon completing the file
    processing the nu_p_sigma_fg field gets populated from the product of
    nu_prompt and sigma_fg. Symbol \f$ \nu_{prompt,g} \f$.
  - NU_DELAYED_BEGIN. Optional. Starts a block that is terminated by a line
    NU_DELAYED_END. Each line in the block processes the first two words as
    [group, nu_delayed]. Populates the nu_delayed field. Upon completing the
    file processing the nu_d_sigma_fg field gets populated from the product of
    nu_delayed and sigma_fg. Symbol \f$ \nu_{delayed,g} \f$.
  - CHI_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_END. Each line in the block processes the first two words as
    [group, chi]. Populates the chi field. Symbol \f$ \chi_{g} \f$.
  - CHI_PROMPT_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_PROMPT_END. Each line in the block processes the first two words as
    [group, chi]. Populates the chi_prompt field.
    Symbol \f$ \chi_{prompt, g} \f$.
  - VELOCITY_BEGIN. Optional. Starts a block that is terminated by a line
    VELOCITY_END. Each line in the block processes the first two words as
    [group, velocity]. Populates the inv_velocity field by inverting parsed
    values. Symbol \f$ \frac{1}{v_g} \f$.
  - INV_VELOCITY_BEGIN. Optional. Starts a block that is terminated by a line
    INV_VELOCITY_END. Each line in the block processes the first two words as
    [group, inv_velocity]. Populates the inv_velocity field. If this field and
    VELOCITY are provided, this field will be used.
    Symbol \f$ \frac{1}{v_g} \f$.
  - PRECURSOR_LAMBDA_BEGIN. Optional. Starts a block that is terminated by a
    line PRECURSOR_LAMBDA_END. Each line in the block processes the first two
    words as [precursor, lambda]. Populates the lambda field (the precursor
    decay constant). Symbol \f$ \lambda_j \f$.
  - PRECURSOR_YIELD_BEGIN. Optional. Starts a block that is terminated by a
    line PRECURSOR_YIELD_END. Each line in the block processes the first two
    words as [precursor, gamma]. Populates the gamma field (the precursor
    production fraction per fission). Symbol \f$ \gamma_j \f$.
  - CHI_DELAYED_BEGIN. Optional. Starts a block that is terminated by a line
    CHI_DELAYED_END. Each line in the block processes the first word as the
    group index and the remaining NUM_PRECURSORS words as the the individual
    precursor's associated delayed spectrum (chi). Populates the chi_d field.
    Symbol \f$ \chi_{delayed,jg} \f$.
  - TRANSFER_MOMENTS_BEGIN. Optional. Starts a block that is terminated by a
    line TRANSFER_MOMENTS_END. Each line in the block processes a line only if
    it starts with the keyword M_GPRIME_G_VAL which needs to be followed by four
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
  MakeFromChiXSFile(const std::string &file_name)
{
  //clear any previous data
  Reset();

  //============================================================
  // Open Chi XS file
  //============================================================

  chi::log.Log()
      << "Reading Chi cross-section file \"" << file_name << "\"\n";

  //opens and checks if open
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open())
    throw std::runtime_error(
        "Failed to open Chi cross-section file "
        "\"" + file_name + "\" in call to " +
        std::string(__FUNCTION__));

  //============================================================
  // Define utility functions for parsing
  //============================================================

  //##################################################
  /// Lambda for reading group structure data.
  auto ReadGroupStructure =
      [](const std::string& keyword,
         std::vector<std::vector<double>>& destination,
         const unsigned int n,
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.assign(n, std::vector<double>(2, 0.0));

        //book-keeping
        std::string line;
        int group;
        double high;
        double low;
        unsigned int count = 0;

        //read the block
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;
        while (line != keyword + "_END")
        {
          //get data from current line
          line_stream >> group >> high >> low;
          destination.at(group).at(0) = high;
          destination.at(group).at(1) = low;
          if (count++ >= n)
            throw std::runtime_error(
                "Too many entries encountered when parsing "
                "group structure.\nThe expected number of entries "
                "is " + std::to_string(n) + ".");

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //##################################################
  /// Lambda for reading vector data.
  auto Read1DData =
      [](const std::string& keyword,
         std::vector<double>& destination,
         const unsigned int n,
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.assign(n, 0.0);

        //book-keeping
        std::string line;
        int i;
        double value;
        unsigned int count = 0;

        //read the bloc
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;
        while (line != keyword + "_END")
        {
          //get data from current line
          line_stream >> i >> value;
          destination.at(i) = value;
          if (count++ >= n)
            throw std::runtime_error(
                "To many entries encountered when parsing "
                "1D data.\nThe expected number of entries is " +
                std::to_string(n));

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //##################################################
  /// Lambda for reading transfer matrix data.
  auto ReadTransferMatrices =
      [](const std::string& keyword,
         std::vector<TransferMatrix>& destination,
         const unsigned int n,  //# of moments
         const unsigned int m,  //# of groups
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.clear();
        for (unsigned int i = 0; i < n; ++i)
          destination.emplace_back(m, m);

        //book-keeping
        std::string word, line;
        double value;
        unsigned int ell;
        unsigned int group;
        unsigned int gprime;

        //read the block
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;
        while (line != keyword + "_END")
        {
          //check that this line contains an entry
          line_stream >> word;
          if (word == "M_GPRIME_G_VAL")
          {
            //get data from current line
            line_stream >> ell >> gprime >> group >> value;
            destination.at(ell).Insert(group, gprime, value);
          }

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //##################################################
  /// Lambda for reading emission spectra data.
  auto ReadEmissionSpectra =
      [](const std::string& keyword,
         EmissionSpectra& destination,
         const unsigned int m, //# of groups
         const unsigned int n, //# of precursors
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.clear();
        for (unsigned int i = 0; i < m; ++i)
          destination.emplace_back(n, 0.0);

        //book-keeping
        std::string word, line;
        double value;
        unsigned int group;
        unsigned int precursor;

        //read the block
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;
        while (line != keyword + "_END")
        {
          //check that this line contains an entry
          line_stream >> word;
          if (word == "G_PRECURSORJ_VAL")
          {
            //get data from current line
            line_stream >> group >> precursor >> value;
            destination.at(group).at(precursor) = value;
          }

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //============================================================
  // Read the Chi XS file
  //============================================================

  std::string word, line;
  unsigned int line_number = 0;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    //parse number of groups
    if (word == "NUM_GROUPS")
    {
      int G;
      line_stream >> G;
      if (G <= 0)
        throw std::logic_error(
            "The specified number of energy groups "
            "must be positive.");
      num_groups = G;
    }

    //parse the number of scattering moments
    if (word == "NUM_MOMENTS")
    {
      int M;
      line_stream >> M;
      if (M < 0)
        throw std::logic_error(
            "The specified number of scattering moments "
            "must be non-negative.");
      scattering_order = std::min(0, M - 1);
    }

    //parse the number of precursors species
    if (word == "NUM_PRECURSORS")
    {
      int J;
      line_stream >> J;
      if (J < 0)
        throw std::logic_error(
            "The specified number of delayed neutron "
            "precursors must be non-negative.");
      num_precursors = J;
    }

    //parse nuclear data
    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f  = file;
      auto& fw = word;

      //============================================================
      // Group Structure Data
      //============================================================

      if (fw == "GROUP_STRUCTURE_BEGIN")
        ReadGroupStructure("GROUP_STRUCTURE",
                           e_bounds, num_groups, f, ls, ln);

      if (fw == "INV_VELOCITY_BEGIN")
        Read1DData("INV_VELOCITY", inv_velocity, num_groups, f, ls, ln);
      if (fw == "VELOCITY_BEGIN" &&
          std::all_of(inv_velocity.begin(), inv_velocity.end(),
                      [](double x) { return x == 0.0; }))
      {
        Read1DData("VELOCITY", inv_velocity, num_groups, f, ls, ln);
        for (unsigned int g = 0; g < num_groups; ++g)
          inv_velocity[g] = 1.0 / inv_velocity[g];
      }

      //==================================================
      // Cross-Section Data
      //==================================================

      if (fw == "SIGMA_T_BEGIN")
        Read1DData("SIGMA_T", sigma_t, num_groups, f, ls, ln);
      if (fw == "SIGMA_A_BEGIN")
        Read1DData("SIGMA_A", sigma_a, num_groups, f, ls, ln);
      if (fw == "SIGMA_F_BEGIN")
        Read1DData("SIGMA_F", sigma_f, num_groups, f, ls, ln);
      if (fw == "NU_SIGMA_F_BEGIN")
        Read1DData("NU_SIGMA_F", nu_sigma_f, num_groups, f, ls, ln);

      //==================================================
      // Neutrons Per Fission
      //==================================================

      if (fw == "NU_BEGIN")
        Read1DData("NU", nu, num_groups, f, ls, ln);
      if (fw == "NU_PROMPT_BEGIN")
        Read1DData("NU_PROMPT", nu_prompt, num_groups, f, ls, ln);
      if (fw == "NU_DELAYED_BEGIN")
        Read1DData("NU_DELAYED", nu_delayed, num_groups, f, ls, ln);
      if (fw == "BETA_BEGIN")
        Read1DData("BETA", beta, num_groups, f, ls, ln);

      //==================================================
      // Fission/Emission Spectra
      //==================================================

      if (fw == "CHI_BEGIN")
        Read1DData("CHI", chi, num_groups, f, ls, ln);
      if (fw == "CHI_PROMPT_BEGIN")
        Read1DData("CHI_PROMPT", chi_prompt, num_groups, f, ls, ln);
      if (num_precursors > 0 && fw == "CHI_DELAYED_BEGIN")
        ReadEmissionSpectra("CHI_DELAYED", chi_delayed,
                            num_groups, num_precursors, f, ls, ln);

      //==================================================
      // Delayed Neutron Precursor Data
      //==================================================

      if (num_precursors > 0)
      {
        if (fw == "PRECURSOR_LAMBDA_BEGIN")
          Read1DData("PRECURSOR_LAMBDA", precursor_lambda,
                     num_precursors, f, ls, ln);

        if (fw == "PRECURSOR_YIELD_BEGIN")
          Read1DData("PRECURSOR_YIELD", precursor_yield,
                     num_precursors, f, ls, ln);
      }

      //==================================================
      // Transfer Data
      //==================================================

      if (fw == "TRANSFER_MOMENTS_BEGIN")
        ReadTransferMatrices("TRANSFER_MOMENTS", transfer_matrices,
                             scattering_order + 1, num_groups,
                             f, ls, ln);

    }//try

    catch (const std::runtime_error& err)
    {
      throw std::runtime_error(
          "Error reading Chi cross-section file "
          "\"" + file_name + "\".\n" +
          "Line number " + std::to_string(line_number) +
          "\n" + err.what());
    }

    catch (const std::out_of_range& err)
    {
      throw std::out_of_range(
          "Error reading Chi cross-section file "
          "\"" + file_name + "\".\n" +
          "Line number " + std::to_string(line_number) +
          "\n" + err.what());
    }

    catch (...)
    {
      throw std::runtime_error(
          "Unknown error encountered in " +
          std::string (__FUNCTION__ ));
    }

    word = "";
  }//while not EOF, read each lines

  //perform checks and enforce physical relationships
  Finalize();


  file.close();
}