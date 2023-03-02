#include "multigroup_xs.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <iostream>

//###################################################################
/**Exports the cross section information to ChiTech format.*/
void chi_physics::MultiGroupXS::
ExportToChiXSFile(const std::string &file_name)
{
  chi::log.Log() << "Exporting transport cross section to file: " << file_name;

  //============================================================
  // Define utility functions
  //============================================================

  /**Lambda to print a 1D-xs*/
  auto Print1DXS = [](std::ofstream& ofile,
                      const std::string& prefix,
                      const std::vector<double>& xs,
                      double min_value=-1.0)
  {
    bool proceed = false;
    if (min_value >= 0.0)
    {
      for (auto val : xs)
        if (val > min_value)
        {
          proceed = true;
          break;
        }

      if (not proceed)
        return;
    }

    ofile << "\n";
    ofile << prefix << "_BEGIN\n";
    {
      unsigned int g = 0;
      for (auto val : xs)
        ofile << g++ << " "
              << val << "\n";
    }
    ofile << prefix << "_END\n";
  };

  //============================================================
  // Open the output file
  //============================================================

  std::ofstream ofile(file_name);

  //============================================================
  // Write the header info
  //============================================================

  std::vector<double> nu, nu_prompt, nu_delayed;
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    if (num_precursors_ > 0)
    {
      nu_prompt.push_back(nu_prompt_sigma_f_[g] / sigma_f_[g]);
      nu_delayed.push_back(nu_delayed_sigma_f_[g] / sigma_f_[g]);
    }
    else
      nu.push_back(nu_sigma_f_[g] / sigma_f_[g]);
  }

  std::vector<double> decay_constants, fractional_yields;
  for (const auto& precursor : precursors_)
  {
    decay_constants.push_back(precursor.decay_constant);
    fractional_yields.push_back(precursor.fractional_yield);
  }



  ofile << "# Exported cross section from ChiTech\n";
  ofile << "# Date: " << chi_objects::ChiTimer::GetLocalDateTimeString() << "\n";
  ofile << "NUM_GROUPS " << num_groups_ << "\n";
  ofile << "NUM_MOMENTS " << scattering_order_ + 1 << "\n";
  if (num_precursors_ > 0)
    ofile << "NUM_PRECURSORS " << num_precursors_ << "\n";

  //basic cross section data
  Print1DXS(ofile, "SIGMA_T", sigma_t_, 1.0e-20);
  Print1DXS(ofile, "SIGMA_A", sigma_a_, 1.0e-20);

  //fission data
  if (!sigma_f_.empty())
  {
    Print1DXS(ofile, "SIGMA_F", sigma_f_, 1.0e-20);
    if (num_precursors_ > 0)
    {
      Print1DXS(ofile, "NU_PROMPT", nu_prompt, 1.0e-20);
      Print1DXS(ofile, "NU_DELAYED", nu_delayed, 1.0e-20);
//      Print1DXS(ofile, "CHI_PROMPT", chi_prompt, 1.0e-20);

      ofile << "\nCHI_DELAYED_BEGIN\n";
      for (unsigned int j = 0; j < num_precursors_; ++j)
        for (unsigned int g = 0; g < num_groups_; ++g)
          ofile << "G_PRECURSOR_VAL"
                << " " << g
                << " " << j
                << " " << precursors_[j].emission_spectrum[g]
                << "\n";
      ofile << "CHI_DELAYED_END\n";

      Print1DXS(ofile, "PRECURSOR_DECAY_CONSTANTS",
                decay_constants, 1.0e-20);
      Print1DXS(ofile, "PRECURSOR_FRACTIONAL_YIELDS",
                fractional_yields, 1.0e-20);

    }
    else
    {
      Print1DXS(ofile, "NU", nu, 1.0e-20);
//      Print1DXS(ofile, "CHI", chi, 1.0e-20);
    }
  }

  //inverse speed data
  if (!inv_velocity_.empty())
    Print1DXS(ofile, "INV_VELOCITY", inv_velocity_, 1.0e-20);

  //transfer matrices
  if (!transfer_matrices_.empty())
  {
    ofile << "\n";
    ofile << "TRANSFER_MOMENTS_BEGIN\n";
    for (size_t ell=0; ell < transfer_matrices_.size(); ++ell)
    {
      if (ell==0) ofile << "#Zeroth moment (l=0)\n";
      else        ofile << "#(l=" << ell << ")\n";

      const auto& matrix = transfer_matrices_[ell];

      for (size_t g=0; g<matrix.rowI_values_.size(); ++g)
      {
        const auto& col_indices = matrix.rowI_indices_[g];
        const auto& col_values  = matrix.rowI_values_[g];

        for (size_t k=0; k<col_indices.size(); ++k)
          ofile << "M_GPRIME_G_VAL "
                << ell << " "
                << col_indices[k] << " "
                << g << " "
                << col_values[k] << "\n";
      }//for g

      ofile << "\n";
    }//for ell
    ofile << "TRANSFER_MOMENTS_END\n";
  }//if has transfer matrices

  if (!production_matrix_.empty())
  {
    ofile << "\n";
    ofile << "PRODUCTION_MATRIX_BEGIN\n";
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      const auto& prod = production_matrix_[g];
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
        ofile << "G_GPRIME_VAL "
              << g << " "
              << gp << " "
              << prod[gp] << "\n";
    }
  }

  ofile.close();

  chi::log.Log0Verbose1() << "Done exporting transport "
                             "cross section to file: " << file_name;
}