#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <string>
#include <algorithm>
#include <numeric>


/**\defgroup ChiXSFile Chi-Tech Cross-section format 1
 *\ingroup LuaPhysics
 *
 * An example Chi-Tech cross section file is shown below. The bare-bones
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
  cross section. Symbol \f$ G \f$.
- NUM_MOMENTS num_moments Required. The number of transfer matrices to allocate
  for this cross section (whether used or not). Typically this number is one
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
  - PRECURSOR_DECAY_CONSTANTS_BEGIN. Optional. Starts a block that is terminated by a
    line PRECURSOR_DECAY_CONSTANTS_END. Each line in the block processes the first two
    words as [precursor, lambda]. Populates the lambda field (the precursor
    decay constant). Symbol \f$ \lambda_j \f$.
  - PRECURSOR_FRACTIONAL_YIELDS_BEGIN. Optional. Starts a block that is terminated by a
    line PRECURSOR_FRACTIONAL_YIELDS_END. Each line in the block processes the first two
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

PRECURSOR_DECAY_CONSTANTS_BEGIN
0		0.1
1   0.2
2   0.3
PRECURSOR_DECAY_CONSTANTS_END

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
G_PRECURSOR_VAL 0  0	1.0
G_PRECURSOR_VAL 0  1	1.0
G_PRECURSOR_VAL 0  2	1.0

G_PRECURSOR_VAL 1  0	0.0
G_PRECURSOR_VAL 1  1	0.0
G_PRECURSOR_VAL 1  2	0.0
CHI_DELAYED_END
\endcode
 * */

//###################################################################
/**This method populates a transport cross section from
 * a Chi cross section file.*/
void chi_physics::TransportCrossSections::
  MakeFromChiXSFile(const std::string &file_name)
{
  Reset();

  //============================================================
  // Open Chi XS file
  //============================================================

  chi::log.Log()
      << "Reading Chi cross section file \"" << file_name << "\"\n";

  //opens and checks if open
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open())
    throw std::runtime_error(
        "Failed to open Chi cross section file "
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
         const unsigned int G,  //# of groups
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.assign(G, std::vector<double>(2, 0.0));

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
          if (count++ >= G)
            throw std::runtime_error(
                "Too many entries encountered when parsing "
                "group structure.\nThe expected number of entries "
                "is " + std::to_string(G) + ".");

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
         const unsigned int N,
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.assign(N, 0.0);

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
          if (count++ >= N)
            throw std::runtime_error(
                "To many entries encountered when parsing "
                "1D data.\nThe expected number of entries is " +
                std::to_string(N));

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
         const unsigned int M,  //# of moments
         const unsigned int G,  //# of groups
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.clear();
        for (unsigned int i = 0; i < M; ++i)
          destination.emplace_back(G, G);

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
         const unsigned int J, //# of precursors
         const unsigned int G, //# of groups
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.clear();
        for (unsigned int j = 0; j < J; ++j)
          destination.emplace_back(G, 0.0);

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
          if (word == "G_PRECURSOR_VAL")
          {
            //get data from current line
            line_stream >> group >> precursor >> value;
            destination.at(precursor).at(group) = value;
          }

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //##################################################
  /// Lambda for reading emission spectra data.
  auto ReadProductionMatrix =
      [](const std::string& keyword,
         EmissionSpectra& destination,
         const unsigned int G, //# of groups
         std::ifstream& file,
         std::istringstream& line_stream,
         unsigned int& line_number)
      {
        //init storage
        destination.clear();
        for (unsigned int g = 0; g < G; ++g)
          destination.emplace_back(G, 0.0);

        //book-keeping
        std::string word, line;
        double value;
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
          if (word == "G_GPRIME_VAL")
          {
            //get data from current line
            line_stream >> group >> gprime >> value;
            destination.at(group).at(gprime) = value;
          }

          //go to next line
          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //##################################################
  /// Lambda for checking for all non-negative values.
  auto is_nonnegative =
      [](const std::vector<double>& vec)
      {
        return !vec.empty() &&
               std::all_of(vec.begin(), vec.end(),
                           [](double x) { return x >= 0.0; });
      };

  //##################################################
  /// Lambda for checking for all strictly positive values.
  auto is_positive =
      [](const std::vector<double>& vec)
      {
        return !vec.empty() &&
               std::all_of(vec.begin(), vec.end(),
                           [](double x) { return x > 0.0; });
      };

  //##################################################
  /// Lambda for checking for any non-zero values.
  auto has_nonzero =
      [](const std::vector<double> & vec)
      {
        return !vec.empty() &&
               std::any_of(vec.begin(), vec.end(),
                           [](double x) { return x > 0.0; });
      };

  //============================================================
  // Read the Chi XS file
  //============================================================

  std::vector<double> decay_constants;
  std::vector<double> fractional_yields;
  std::vector<std::vector<double>> emission_spectra;
  std::vector<double> nu, nu_prompt, nu_delayed, beta;

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
            "must be strictly positive.");
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
      scattering_order = std::max(0, M - 1);
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
      precursors.resize(num_precursors);
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
      {
        Read1DData("INV_VELOCITY", inv_velocity, num_groups, f, ls, ln);
        if (!is_positive(inv_velocity))
          throw std::logic_error(
              "Invalid inverse velocity value encountered.\n"
              "Only strictly positive values are permitted.");
      }
      if (fw == "VELOCITY_BEGIN" && inv_velocity.empty())
      {
        Read1DData("VELOCITY", inv_velocity, num_groups, f, ls, ln);
        if (!is_positive(inv_velocity))
          throw std::logic_error(
              "Invalid velocity value encountered.\n"
              "Only strictly positive values are permitted.");

        //compute inverse
        for (unsigned int g = 0; g < num_groups; ++g)
          inv_velocity[g] = 1.0 / inv_velocity[g];
      }

      //==================================================
      // Cross Section Data
      //==================================================

      if (fw == "SIGMA_T_BEGIN")
      {
        Read1DData("SIGMA_T", sigma_t, num_groups, f, ls, ln);
        if (!is_nonnegative(sigma_t))
          throw std::logic_error(
              "Invalid total cross section value encountered.\n"
              "Negative values are not permitted.");
      }//if sigma_t

      if (fw == "SIGMA_A_BEGIN")
      {
        Read1DData("SIGMA_A", sigma_a, num_groups, f, ls, ln);
        if (!is_nonnegative(sigma_a))
          throw std::logic_error(
              "Invalid absorption cross section value encountered.\n"
              "Negative values are not permitted.");
      }//if sigma_a

      if (fw == "SIGMA_F_BEGIN")
      {
        Read1DData("SIGMA_F", sigma_f, num_groups, f, ls, ln);
        if (!has_nonzero(sigma_f))
        {
          chi::log.Log0Warning()
              << "The fission cross section specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          sigma_f.clear();
        }
        if (!is_nonnegative(sigma_f))
          throw std::logic_error(
              "Invalid fission cross section value encountered.\n"
              "Negative values are not permitted.");
      }//if sigma_f

      if (fw == "NU_SIGMA_F_BEGIN")
      {
        Read1DData("NU_SIGMA_F", nu_sigma_f, num_groups, f, ls, ln);
        if (!has_nonzero(nu_sigma_f))
        {
          chi::log.Log0Warning()
              << "The production cross-section specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          nu_sigma_f.clear();
        }
        if (!is_nonnegative(nu_sigma_f))
          throw std::logic_error(
              "Invalid production cross section value encountered.\n"
              "Negative values are not permitted.");
      }//if nu_sigma_f

      //==================================================
      // Neutrons Per Fission
      //==================================================

      if (fw == "NU_BEGIN")
      {
        Read1DData("NU", nu, num_groups, f, ls, ln);
        if (!has_nonzero(nu))
        {
          chi::log.Log0Warning()
              << "The total fission neutron yield specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          nu.clear();
        }
        if (!std::all_of(nu.begin(), nu.end(),
                         [](double x)
                         { return x == 0.0 || x > 1.0; }))
          throw std::logic_error(
              "Invalid total fission neutron yield value encountered.\n"
              "Only values strictly greater than one, or zero, are "
              "permitted.");

        //compute prompt/delayed nu, if needed
        if (num_precursors > 0 &&
            !nu.empty() && !beta.empty() &&
            nu_prompt.empty() && nu_delayed.empty())
        {
          nu_prompt.assign(num_groups, 0.0);
          nu_delayed.assign(num_groups, 0.0);
          for (unsigned int g = 0; g < num_groups; ++g)
          {
            nu_prompt[g] = (1.0 - beta[g]) * nu[g];
            nu_delayed[g] = beta[g] * nu[g];
          }
        }
      }//if nu

      if (fw == "NU_PROMPT_BEGIN")
      {
        Read1DData("NU_PROMPT", nu_prompt, num_groups, f, ls, ln);
        if (!has_nonzero(nu_prompt))
        {
          chi::log.Log0Warning()
              << "The prompt fission neutron yield specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          nu_prompt.clear();
        }
        if (!std::all_of(nu_prompt.begin(), nu_prompt.end(),
                         [](double x)
                         { return x == 0.0 || x > 1.0; }))
          throw std::logic_error(
              "Invalid prompt fission neutron yield value encountered.\n"
              "Only values strictly greater than one, or zero, are "
              "permitted.");
      }

      if (fw == "NU_DELAYED_BEGIN")
      {
        Read1DData("NU_DELAYED", nu_delayed, num_groups, f, ls, ln);
        if (!has_nonzero(nu_delayed))
        {
          chi::log.Log0Warning()
              << "The delayed fission neutron yield specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          nu_prompt.clear();
        }
        if (!is_nonnegative(nu_delayed))
          throw std::logic_error(
              "Invalid delayed fission neutron yield value encountered.\n"
              "Only non-negative values are permitted.");
      }

      if (fw == "BETA_BEGIN")
      {
        Read1DData("BETA", beta, num_groups, f, ls, ln);
        if (!has_nonzero(beta))
        {
          chi::log.Log0Warning()
              << "The delayed neutron fraction specified in "
              << "\"" << file_name << "\" is uniformly zero..."
              << "Clearing it.";
          beta.clear();
        }
        if (!std::all_of(beta.begin(), beta.end(),
                         [](double x)
                         { return x >= 0.0 && x <= 1.0; }))
          throw std::logic_error(
              "Invalid delayed neutron fraction value encountered.\n"
              "Only values in the range [0.0, 1.0] are permitted.");

        //compute prompt/delayed nu, if needed
        if (num_precursors > 0 &&
            !nu.empty() && !beta.empty() &&
            nu_prompt.empty() && nu_delayed.empty())
        {
          nu_prompt.assign(num_groups, 0.0);
          nu_delayed.assign(num_groups, 0.0);
          for (unsigned int g = 0; g < num_groups; ++g)
          {
            nu_prompt[g] = (1.0 - beta[g]) * nu[g];
            nu_delayed[g] = beta[g] * nu[g];
          }
        }
      }//if beta

      //==================================================
      // Fission/Emission Spectra
      //==================================================

      if (fw == "CHI_BEGIN")
      {
        Read1DData("CHI", chi, num_groups, f, ls, ln);
        if (!has_nonzero(chi))
          throw std::logic_error(
              "Invalid steady-state fission spectrum encountered.\n"
              "No non-zero values found.");
        if (!is_nonnegative(chi))
          throw std::logic_error(
              "Invalid steady-state fission spectrum value encountered.\n"
              "Only non-negative values are permitted.");

        //normalizing
        double sum = std::accumulate(chi.begin(), chi.end(), 0.0);
        std::transform(chi.begin(), chi.end(), chi.begin(),
                      [sum](double& x) { return x / sum; });
      }//if chi

      if (fw == "CHI_PROMPT_BEGIN")
      {
        Read1DData("CHI_PROMPT", chi_prompt, num_groups, f, ls, ln);
        if (!has_nonzero(chi_prompt))
          throw std::logic_error(
              "Invalid prompt fission spectrum encountered.\n"
              "No non-zero values found.");
        if (!is_nonnegative(chi_prompt))
          throw std::logic_error(
              "Invalid prompt fission spectrum value encountered.\n"
              "Only non-negative values are permitted.");

        //normalizing
        double sum = std::accumulate(chi_prompt.begin(), chi_prompt.end(), 0.0);
        std::transform(chi_prompt.begin(),
                       chi_prompt.end(),
                       chi_prompt.begin(),
                       [sum](double& x) { return x / sum; });

      }//if prompt chi

      if (num_precursors > 0 &&
          fw == "CHI_DELAYED_BEGIN")
      {
        ReadEmissionSpectra("CHI_DELAYED", emission_spectra,
                            num_precursors, num_groups, f, ls, ln);
        for (unsigned int j = 0; j < num_precursors; ++j)
        {
          if (!has_nonzero(emission_spectra[j]))
            throw std::logic_error(
                "Invalid delayed emission spectrum encountered for "
                "precursor species " + std::to_string(j) + ".\n" +
                "No non-zero values found.");
          if (!is_nonnegative(emission_spectra[j]))
            throw std::logic_error(
                "Invalid delayed emission spectrum value encountered "
                "for precursor species " + std::to_string(j) + ".\n" +
                "Only non-negative values are permitted.");

          //normalizing
          double sum = std::accumulate(emission_spectra[j].begin(),
                                       emission_spectra[j].end(), 0.0);
          std::transform(emission_spectra[j].begin(),
                         emission_spectra[j].end(),
                         emission_spectra[j].begin(),
                         [sum](double& x) { return x / sum; });
        }
      }//if delayed chi

      //==================================================
      // Delayed Neutron Precursor Data
      //==================================================

      if (num_precursors > 0)
      {
        if (fw == "PRECURSOR_DECAY_CONSTANTS_BEGIN")
        {
          Read1DData("PRECURSOR_DECAY_CONSTANTS",
                     decay_constants, num_precursors, f, ls, ln);
          if (!is_positive(decay_constants))
            throw std::logic_error(
                "Invalid precursor decay constant value encountered.\n"
                "Only strictly positive values are permitted.");
        }//if decay constants

        if (fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN")
        {
          Read1DData("PRECURSOR_FRACTIONAL_YIELDS",
                     fractional_yields, num_precursors, f, ls, ln);
          if (!has_nonzero(fractional_yields))
            throw std::logic_error(
                "Invalid precursor fractional yields encountered.\n"
                "No non-zero values found.");
          if (!std::all_of(fractional_yields.begin(),
                           fractional_yields.end(),
                           [](double x) { return x >= 0.0 && x <= 1.0; }))
            throw std::logic_error(
                "Invalid precursor fractional yield value encountered.\n"
                "Only values in the range [0.0, 1.0] are permitted.");

          //normalizing
          double sum = std::accumulate(fractional_yields.begin(),
                                       fractional_yields.end(), 0.0);
          std::transform(fractional_yields.begin(),
                         fractional_yields.end(),
                         fractional_yields.begin(),
                         [sum](double& x) { return x / sum; });


        }
      }

      //==================================================
      // Transfer Data
      //==================================================

      if (fw == "TRANSFER_MOMENTS_BEGIN")
        ReadTransferMatrices("TRANSFER_MOMENTS", transfer_matrices,
                             scattering_order + 1, num_groups, f, ls, ln);

      if (fw == "PRODUCTION_MATRIX_BEGIN")
        ReadProductionMatrix("PRODUCTION_MATRIX",
                             production_matrix, num_groups, f, ls, ln);
    }//try

    catch (const std::runtime_error& err)
    {
      throw std::runtime_error(
          "Error reading Chi cross section file "
          "\"" + file_name + "\".\n" +
          "Line number " + std::to_string(line_number) +
          "\n" + err.what());
    }

    catch (const std::logic_error& err)
    {
      throw std::logic_error(
          "Error reading Chi cross section file "
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
  file.close();

  //============================================================
  // Compute auxiliary data
  //============================================================

  if (sigma_a.empty())
    ComputeAbsorption();

  ComputeDiffusionParameters();

  //============================================================
  // Compute and check fission data
  //============================================================

  //determine if the material is fissionable
  is_fissionable = !sigma_f.empty() || !nu_sigma_f.empty();

  //clear fission data if not fissionable
  if (!is_fissionable)
  {
    sigma_f.clear();
    nu_sigma_f.clear();
    nu_prompt_sigma_f.clear();
    nu_delayed_sigma_f.clear();
    chi.clear();
    chi_prompt.clear();
    precursors.clear();
  }//if not fissionable

  //otherwise, set fission data based on inputs
  else
  {
    //set prompt/delayed fission data
    if (num_precursors > 0)
    {
      //check spectra
      if (chi_prompt.empty() ||
          std::any_of(emission_spectra.begin(),
                      emission_spectra.end(),
                      [](const std::vector<double>& x)
                      { return x.empty(); }))
        throw std::logic_error(
            "Invalid fission/emission spectrum encountered.\n"
            "Either the prompt fission spectrum or a delayed "
            "emission spectrum was not found.");

      //check fission neutron yields
      if (nu_prompt.empty() && nu_delayed.empty())
        throw std::logic_error(
            "Invalid fission data specification encountered.\n"
            "Either prompt and delay fission neutron yields or "
            "the total fission neutron yield and delayed neutron "
            "fraction must be specified.");

      //compute the fission cross section
      if (sigma_f.empty())
      {
        sigma_f.assign(num_groups, 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          if (nu_sigma_f[g] > 0.0)
          {
            double nu_total = nu_prompt[g] + nu_delayed[g];
            sigma_f[g] = nu_sigma_f[g] / nu_total;
          }
      }

      //compute production cross sections
      nu_sigma_f.assign(num_groups, 0.0);
      nu_prompt_sigma_f.assign(num_groups, 0.0);
      nu_delayed_sigma_f.assign(num_groups, 0.0);
      for (unsigned int g = 0; g < num_groups; ++g)
      {
        double nu_total = nu_prompt[g] + nu_delayed[g];
        nu_sigma_f[g] = nu_total * sigma_f[g];
        nu_prompt_sigma_f[g] = nu_prompt[g] * sigma_f[g];
        nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
      }

      //add data to the precursor structs
      for (unsigned int j = 0; j < num_precursors; ++j)
      {
        precursors[j].decay_constant = decay_constants[j];
        precursors[j].fractional_yield = fractional_yields[j];
        precursors[j].emission_spectrum = emission_spectra[j];
      }

      //create production matrix, if empty
      if (production_matrix.empty())
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          std::vector<double> vals;
          for (unsigned int gp = 0; gp < num_groups; ++gp)
            vals.push_back(chi_prompt[g] * nu_prompt_sigma_f[gp]);
          production_matrix.push_back(vals);
        }
    }//prompt/delayed

    //set steady-state fission data
    else
    {
      //check spectra
      if (chi.empty())
        throw std::logic_error(
            "Invalid steady-state fission spectrum encountered.\n"
            "No steady-state fission spectrum found.");

      //check nu
      if (nu.empty())
        throw std::logic_error(
            "Invalid total fission neutron yield encountered.\n"
            "No total fission neutron yield found.");

      //compute fission/production cross section
      if (!sigma_f.empty())
      {
        sigma_f.assign(num_groups, 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          if (nu_sigma_f[g] > 0.0)
            sigma_f[g] = nu_sigma_f[g] / nu[g];
      }
      else
      {
        nu_sigma_f.assign(num_groups, 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          if (sigma_f[g] > 0.0)
            nu_sigma_f[g] = nu[g] * sigma_f[g];
      }

      //create production matrix, if empty
      if (production_matrix.empty())
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          std::vector<double> vals;
          for (unsigned int gp = 0; gp < num_groups; ++gp)
            vals.push_back(chi[g] * nu_sigma_f[gp]);
          production_matrix.push_back(vals);
        }
    }//steady-state fission
  }//if fissionable
}
