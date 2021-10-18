#include "material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <string>



void chi_physics::TransportCrossSections::FinalizeCrossSections()
{
  const double eps = 1.0e-8;

  //======================================== Estimates sigma_a group-by-group
  double sigma_a_sum = 0.0;
  for (size_t g = 0; g < num_groups; ++g)
    sigma_a_sum += sigma_a[g];

  if (not transfer_matrices.empty() and (sigma_a_sum < 1.0e-28))
  {
    sigma_a = ComputeAbsorptionXSFromTransfer();

    chi_log.Log(LOG_0WARNING)
        << __FUNCTION__
        << ": sigma_a was estimated from the transfer matrix.";
  }

  //======================================== Determine if fissile or not
  double sigma_f_sum = 0.0;
  for (int g = 0; g < num_groups; ++g)
    sigma_f_sum += sigma_f[g];

  if (sigma_f_sum > 1.0e-28) is_fissile = true;
  if (not is_fissile and num_precursors > 0)
  {
    num_precursors = 0;
    precursor_lambda.clear();
    precursor_yield.clear();
    chi_delayed.clear();

    chi_log.Log(LOG_ALLWARNING)
        << __FUNCTION__
        << ": Precursors found in a non-fissile material. "
        << "Setting the number of precursors to zero for consistency.";
  }

  if (is_fissile)
  {
    //======================================== Determine present nu terms
    bool has_nu = false;
    bool has_nu_p = false;
    bool has_nu_d = false;
    for (int g = 0; g < num_groups; ++g)
    {
      if (not has_nu and nu[g] > 0.0) has_nu = true;
      if (not has_nu_p and nu_prompt[g] > 0.0) has_nu_p = true;
      if (not has_nu_d and nu_delayed[g] > 0.0) has_nu_d = true;
    }

    //======================================== Determine present chi terms
    bool has_chi = false;
    bool has_chi_p = false;
    for (int g = 0; g < num_groups; ++g)
    {
      if (not has_chi and chi[g] > 0.0) has_chi = true;
      if (not has_chi_p and chi_prompt[g] > 0.0) has_chi_p = true;
    }

    bool has_chi_d = false;
    if (num_precursors > 0)
    {
      for (int j = 0; j < num_precursors; ++j)
      {
        double chi_dj_sum = 0.0;
        for (int g = 0; g < num_groups; ++g)
          chi_dj_sum += chi_delayed[g][j];

        if (chi_dj_sum < eps)
        {
          chi_log.Log(LOG_ALLERROR)
            << __FUNCTION__ << ": Precursor family " << j
            << "does not have a non-zero fission spectrum.";
          exit(EXIT_FAILURE);
        }
      }
      has_chi_d = true;
    }//if num_precursors > 0

    //======================================== Check precursor terms
    if (num_precursors > 0)
    {
      for (int j = 0; j < num_precursors; ++j)
      {
        if (precursor_lambda[j] < eps)
        {
          chi_log.Log(LOG_ALLERROR)
            << __FUNCTION__ << ": Precursor family " << j << "'s decay "
            << "constant must be non-zero.";
          exit(EXIT_FAILURE);
        }
        if (precursor_yield[j] < eps)
        {
          chi_log.Log(LOG_ALLERROR)
              << __FUNCTION__ << ": Precursor family " << j << "'s yield "
              << "must be non-zero.";
          exit(EXIT_FAILURE);
        }
      }
    }

    //======================================== Normalization
    if (has_chi)
    {
      double chi_sum = 0.0;
      for (auto& v : chi) chi_sum += v;
      if (fabs(chi_sum - 1.0) > eps and chi_sum > eps)
      {
        chi_log.Log(LOG_ALLWARNING)
          << __FUNCTION__ << ": Total fission spectrum does "
          << "not sum to unity. Normalizing.";
        for (auto& v : chi) v /= chi_sum;
      }
    }

    if (has_chi_p)
    {
      double chi_p_sum = 0.0;
      for (auto& v : chi_prompt) chi_p_sum += v;
      if (fabs(chi_p_sum - 1.0) > eps and chi_p_sum > eps)
      {
        chi_log.Log(LOG_ALLWARNING)
            << __FUNCTION__ << ": Prompt fission spectrum does "
            << "not sum to unity. Normalizing.";
        for (auto& v : chi_prompt) v /= chi_p_sum;
      }
    }

    if (has_chi_d)
    {
      for (int j = 0; j < num_precursors; ++j)
      {
        double chi_dj_sum = 0.0;
        for (int g = 0; g < num_groups; ++g)
          chi_dj_sum += chi_delayed[g][j];
        if (fabs(chi_dj_sum - 1.0) > eps and chi_dj_sum > eps)
        {
          chi_log.Log(LOG_ALLWARNING)
              << __FUNCTION__ << ": Delayed fission spectrum for precursor "
              << "family " << j << " does not sum to unity. Normalizing.";
          for (int g = 0; g < num_groups; ++g)
            chi_delayed[g][j] /= chi_dj_sum;
        }
      }
    }

    if (num_precursors > 0)
    {
      double yield_sum = 0.0;
      for (auto& v : precursor_yield) yield_sum += v;
      if (fabs(yield_sum - 1.0) > eps and yield_sum > eps)
      {
        chi_log.Log(LOG_ALLWARNING)
            << __FUNCTION__ << ": Precursor yield sum does not sum  "
            << "to unity. Normalizing.";
        for (auto& v : precursor_yield) v /= yield_sum;
      }
    }

    //======================================== Check for input compatibility
    if ((has_nu_p and not has_chi_p) or (not has_nu_p and has_chi_p))
    {
      chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": If prompt nu or chi are provided, the other "
        << "must be as well.";
      exit(EXIT_FAILURE);
    }
    if ((has_nu_d and not has_chi_d) or (not has_nu_d and has_chi_d))
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": If prompt nu or chi are provided, the other "
          << "must be as well.";
      exit(EXIT_FAILURE);
    }
    if ((has_nu and not has_chi) or (not has_nu and has_chi))
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": If total nu or chi are provided, the other "
          << "must be as well.";
      exit(EXIT_FAILURE);
    }

    //============================== Compute total from prompt and delayed
    if (has_nu_p and has_nu_d)
    {
      for (int g = 0; g < num_groups; ++g)
        nu[g] = nu_prompt[g] + nu_delayed[g];
      has_nu = true;
    }

    if (has_chi_p and has_chi_d)
    {
      for (int g = 0; g < num_groups; ++g)
      {
        double beta = nu_delayed[g] / nu[g];
        chi[g] = (1.0 - beta) * chi_prompt[g];
        for (int j = 0; j < num_precursors; ++j)
          chi[g] += beta * precursor_yield[j] * chi_delayed[g][j];
      }
      has_chi = true;
    }

    //======================================== Make sure appropriate quantities
    //                                         were provided
    if (num_precursors > 0)
    {
      if (not has_nu_p or not has_nu_d)
      {
        chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": Both prompt and delayed nu must be provided "
          << "when precursors are present.";
        exit(EXIT_FAILURE);
      }
      if (not has_chi_p or not has_chi_d)
      {
        chi_log.Log(LOG_ALLERROR)
            << __FUNCTION__ << ": Both prompt and delayed chi must be "
            << "provided when precursors are present.";
        exit(EXIT_FAILURE);
      }
    }
    else
    {
      if (not has_nu)
      {
        chi_log.Log(LOG_ALLERROR)
            << __FUNCTION__ << ": Total nu must be provided "
            << "when precursors are present.";
        exit(EXIT_FAILURE);
      }
      if (not has_chi)
      {
        chi_log.Log(LOG_ALLERROR)
            << __FUNCTION__ << ": Total chi must be "
            << "provided when precursors are present.";
        exit(EXIT_FAILURE);
      }
    }

    //compute nu_sigma_f terms
    for (size_t g = 0; g < num_groups; ++g)
    {
      nu_sigma_f        [g] = nu        [g] * sigma_f[g];
      nu_prompt_sigma_f [g] = nu_prompt [g] * sigma_f[g];
      nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
    }

  }//if fissile
}
