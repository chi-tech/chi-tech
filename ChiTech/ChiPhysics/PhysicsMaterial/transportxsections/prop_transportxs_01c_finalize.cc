#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <numeric>
#include <algorithm>
#include <string>

void chi_physics::TransportCrossSections::FinalizeCrossSections()
{
  const double eps = 1.0e-12;

  //============================================================
  // Set the absorption cross-section, if unset
  //============================================================

  if (std::all_of(sigma_a.begin(), sigma_a.end(),
                  [](double x) { return x == 0.0; }))
    ComputeAbsorption();

  //============================================================
  // Define utility functions
  //============================================================

  auto was_xs_specified =
      [](const std::vector<double>& vec)
      {
        return std::all_of(vec.begin(), vec.end(),
                           [](double x) { return x > 0.0; });
      };

  auto was_spectrum_specified =
      [](const std::vector<double>& vec)
      {
        return std::all_of(vec.begin(), vec.end(),
                           [](double x) { return x >= 0.0; });
      };

  //============================================================
  // Determine if fissionable or not
  //============================================================

  is_fissionable = was_xs_specified(sigma_f) ||
                   was_xs_specified(nu_sigma_f);

  //============================================================
  // Zero fission data if not fissionable
  //============================================================

  if (!is_fissionable)
  {
    num_precursors = 0;

    sigma_f.assign(num_groups, 0.0);
    nu_sigma_f.assign(num_groups, 0.0);
    nu_prompt_sigma_f.assign(num_groups, 0.0);
    nu_delayed_sigma_f.assign(num_groups, 0.0);

    nu.assign(num_groups, 0.0);
    nu_prompt.assign(num_groups, 0.0);
    nu_delayed.assign(num_groups, 0.0);

    chi.assign(num_groups, 0.0);
    chi_prompt.assign(num_groups, 0.0);
    chi_delayed.clear();

    precursor_lambda.clear();
    precursor_yield.clear();
  }

  //============================================================
  // Check specified fission data
  //============================================================

  else
  {
    try
    {
      //==================================================
      // Check prompt/delayed specification
      //==================================================

      if (num_precursors > 0)
      {
        //checks for cross-section specification
        if (was_xs_specified(sigma_f) &&
            was_xs_specified(nu_prompt) &&
            was_xs_specified(nu_delayed))
        {
          if (!std::all_of(nu_prompt.begin(), nu_prompt.end(),
                          [](double x) { return x > 1.0; }))
            throw std::runtime_error(
                "Prompt fission must yield more than one "
                "neutron on average.");

          //compute other properties
          for (unsigned int g = 0; g < num_groups; ++g)
          {
            nu[g] = nu_prompt[g] + nu_delayed[g];
            nu_sigma_f[g] = nu[g] * sigma_f[g];
            nu_prompt_sigma_f[g] = nu_prompt[g] * sigma_f[g];
            nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
          }
        }
        else
          throw std::runtime_error(
              "Invalid specification of fission data.");

        //check and normalize prompt fission spectrum
        {
          //throw error if not specified
          if (!was_spectrum_specified(chi_prompt))
            throw std::runtime_error(
                "The prompt fission spectrum was not provided.");

          //normalize the prompt fission spectrum
          double chip_sum = std::accumulate(chi_prompt.begin(),
                                            chi_prompt.end(), 0.0);
          for (unsigned int g = 0; g < num_groups; ++g)
            chi_prompt[g] /= chip_sum;
        }

        //check and normalize delayed emission spectra
        {
          //create chi_delayed transpose for easier checks
          EmissionSpectra tmp;
          tmp.resize(num_precursors, std::vector<double>(num_groups, 0.0));
          for (unsigned int g = 0; g < num_groups; ++g)
            for (unsigned int j = 0; j < num_precursors; ++j)
              tmp[j][g] = chi_delayed[g][j];

          for (unsigned int j = 0; j < num_precursors; ++j)
          {
            //throw error if emission spectrum j was not specified
            if (!was_spectrum_specified(tmp[j]))
              throw std::runtime_error(
                  "The delayed emission spectrum for precursor "
                  "species " + std::to_string(j) + " was not provided.");

            //normalize each emission spectrum
            double chidj_sum = std::accumulate(tmp[j].begin(),
                                               tmp[j].end(), 0.0);
            for (unsigned int g = 0; g < num_groups; ++g)
              chi_delayed[g][j] /= chidj_sum;
          }//for j
        }

        //check the precursor data
        {
          //throw error if decay constants were not specified
          if (!was_xs_specified(precursor_lambda))
            throw std::runtime_error(
                "The delayed neutron precursor decay constants "
                "was not provided.");

          //throw error if not specified
          if (!was_spectrum_specified(precursor_yield))
            throw std::runtime_error(
                "The delayed neutron precursor yields was not provided.");

          //normalize the precursor yields
          double yield_sum = std::accumulate(precursor_yield.begin(),
                                             precursor_yield.end(), 0.0);
          for (unsigned int j = 0; j < num_precursors; ++j)
            precursor_yield[j] /= yield_sum;
        }

        //compute beta-weighted total fission spectrum
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          //compute beta-averaged total fission spectrum
          double beta = nu_delayed[g] / nu[g];
          chi[g] = (1.0 - beta) * chi_prompt[g];
          for (unsigned int j = 0; j < num_precursors; ++j)
            chi[g] += beta * precursor_yield[j] * chi_delayed[g][j];
        }

        //normalize total chi just in case
        double chi_sum = std::accumulate(chi.begin(), chi.end(), 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          chi[g] /= chi_sum;
      }//if prompt/delayed

      //==================================================
      // Check total fission specification
      //==================================================
      else
      {
        //check that nu was specified correctly
        if (!was_xs_specified(nu))
          throw std::runtime_error(
              "Total neutrons per fission was not provided.");
        if (!std::all_of(nu.begin(), nu.end(),
                         [](double x) { return x > 1.0; }))
          throw std::runtime_error(
              "Total fission yield must be greater than unity.");

        //compute other quantities
        if (was_xs_specified(sigma_f))
          for (unsigned int g = 0; g < num_groups; ++g)
            nu_sigma_f[g] = nu[g] * sigma_f[g];
        else if (was_xs_specified(nu_sigma_f))
          for (unsigned int g = 0; g < num_groups; ++g)
            sigma_f[g] = nu_sigma_f[g] / nu[g];
        else
          throw std::runtime_error(
              "Neither the fission cross-section nor the "
              "fission multiplicity cross-section was specified.");

        //check that chi was specified
        if (!was_spectrum_specified(chi))
          throw std::runtime_error(
              "Total fission spectrum was not provided.");

        //normalize the total fission spectrum
        double chi_sum = std::accumulate(chi.begin(), chi.end(), 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          chi[g] /= chi_sum;
      }
    }
    catch (const std::runtime_error& err)
    {
      chi::log.LogAllError()
          << "Invalid fission data encountered.\n"
          << "Error Message:\n" << err.what();
      chi::Exit(EXIT_FAILURE);
    }
    catch (...)
    {
      chi::log.LogAllError()
          << "Unexpected error encountered in " << __FUNCTION__ ;
      chi::Exit(EXIT_FAILURE);
    }
  }

  //============================================================
  // Compute diffusion parameters
  //============================================================

  ComputeDiffusionParameters();
}
