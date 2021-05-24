#include "material_property_transportxsections.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <iostream>

//###################################################################
/**Exports the cross-section information to ChiTech format.*/
void chi_physics::TransportCrossSections::
  ExportToChiFormat(const std::string &file_name)
{
  ChiLog& chi_log = ChiLog::GetInstance();
  chi_log.Log() << "Exporting transport cross section to file: " << file_name;

  std::ofstream ofile(file_name);

  //======================================== Writing header info
  ofile << "# Exported cross section from ChiTech\n";
  ofile << "# Date: " << ChiTimer::GetLocalDateTimeString() << "\n";
  ofile << "NUM_GROUPS " << num_groups << "\n";
  ofile << "NUM_MOMENTS " << scattering_order+1 << "\n";
  if (num_precursors>0)
    ofile << "NUM_PRECURSORS " << num_precursors << "\n";

  //======================================== Sigma_t
  ofile << "\n";
  ofile << "SIGMA_T_BEGIN\n";
  {
    int g=0;
    for (auto val : sigma_tg)
      ofile << g++ << " " << val << "\n";
  }
  ofile << "SIGMA_T_END\n";

  //======================================== Sigma_f
  bool sigma_f_filled = false;
  for (auto val : sigma_fg)
    if (val > 1.0e-20) sigma_f_filled = true;

  if (sigma_f_filled)
  {
    ofile << "\n";
    ofile << "SIGMA_F_BEGIN\n";
    {
      int g=0;
      for (auto val : sigma_tg)
        ofile << g++ << " " << val << "\n";
    }
    ofile << "SIGMA_F_END\n";
  }

  //======================================== Nu
  bool nu_filled = false;
  for (auto val : nu_sigma_fg)
    if (val > 1.0e-20) nu_filled = true;

  if (nu_filled)
  {
    ofile << "\n";
    ofile << "NU_BEGIN\n";
    {
      int g=0;
      for (auto val : nu_sigma_fg)
      {
        int gval = g++;
        ofile << gval << " " << val/sigma_fg[g-1] << "\n";
      }
    }
    ofile << "NU_END\n";
  }

  //======================================== Nu-prompt
  bool nu_prompt_filled = false;
  for (auto val : nu_p_sigma_fg)
    if (val > 1.0e-20) nu_prompt_filled = true;

  if (nu_prompt_filled)
  {
    ofile << "\n";
    ofile << "NU_PROMPT_BEGIN\n";
    {
      int g=0;
      for (auto val : nu_p_sigma_fg)
      {
        int gval = g++;
        ofile << gval << " " << val/sigma_fg[g-1] << "\n";
      }
    }
    ofile << "NU_PROMPT_END\n";
  }

  //======================================== Nu-delayed
  bool nu_delayed_filled = false;
  for (auto val : nu_d_sigma_fg)
    if (val > 1.0e-20) nu_delayed_filled = true;

  if (nu_delayed_filled)
  {
    ofile << "\n";
    ofile << "NU_DELAYED_BEGIN\n";
    {
      int g=0;
      for (auto val : nu_d_sigma_fg)
      {
        int gval = g++;
        ofile << gval << " " << val/sigma_fg[g-1] << "\n";
      }
    }
    ofile << "NU_DELAYED_END\n";
  }

  //======================================== Chi-prompt
  bool chi_prompt_filled = false;
  for (auto val : chi_g)
    if (val > 1.0e-20) chi_prompt_filled = true;

  if (chi_prompt_filled)
  {
    ofile << "\n";
    ofile << "CHI_PROMPT_BEGIN\n";
    {
      int g=0;
      for (auto val : chi_g)
      {
        int gval = g++;
        ofile << gval << " " << val << "\n";
      }
    }
    ofile << "CHI_PROMPT_END\n";
  }

  //======================================== Chi-delayed

  //======================================== ddt-coeff
  bool ddt_filled = false;
  for (auto val : ddt_coeff)
    if (val > 1.0e-20) ddt_filled = true;

  if (ddt_filled)
  {
    ofile << "\n";
    ofile << "DDT_COEFF_BEGIN\n";
    {
      int g=0;
      for (auto val : ddt_coeff)
        ofile << g++ << " " << val << "\n";
    }
    ofile << "DDT_COEFF_END\n";
  }

  //======================================== Transfer matrices
  if (not transfer_matrix.empty())
  {
    ofile << "TRANSFER_MOMENTS_BEGIN\n";
    for (int ell=0; ell<transfer_matrix.size(); ++ell)
    {
      if (ell==0) ofile << "#Zeroth moment (l=0)\n";
      else        ofile << "#(l=" << ell << ")\n";

      const auto& matrix = transfer_matrix[ell];

      for (int g=0; g<matrix.rowI_values.size(); ++g)
      {
        const auto& col_indices = matrix.rowI_indices[g];
        const auto& col_values  = matrix.rowI_values[g];

        for (int k=0; k<col_indices.size(); ++k)
        {
          size_t gprime = col_indices[k];
          ofile << "M_GPRIME_G_VAL " << ell << " "
                                     << gprime << " "
                                     << g << " "
                                     << col_values[k] << "\n";
        }
      }//for g

      ofile << "\n";
    }//for ell
    ofile << "TRANSFER_MOMENTS_END\n";
  }//if has transfer matrices

  ofile.close();

  chi_log.Log(LOG_0VERBOSE_1) << "Done exporting transport "
                                 "cross section to file: " << file_name;
}