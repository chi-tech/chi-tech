#include "multigroup_xs.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>


//######################################################################
void chi_physics::MultiGroupXS::Reset()
{
  num_groups_ = 0;
  scattering_order_ = 0;
  num_precursors_ = 0;
  is_fissionable_ = false;

  sigma_t_.clear();
  sigma_f_.clear();
  sigma_a_.clear();

  nu_sigma_f_.clear();
  nu_prompt_sigma_f_.clear();
  nu_delayed_sigma_f_.clear();

  inv_velocity_.clear();

  transfer_matrices_.clear();
  production_matrix_.clear();

  precursors_.clear();

  //Diffusion quantities
  diffusion_initialized_ = false;
  diffusion_coeff_.clear();
  sigma_removal_.clear();
  sigma_s_gtog_.clear();

  //Monte-Carlo quantities
  scattering_initialized_ = false;
  cdf_gprime_g_.clear();
  scat_angles_gprime_g_.clear();
}


//######################################################################
/** Makes a simple material with no transfer matrix just sigma_t. */
void chi_physics::MultiGroupXS::
MakeSimple0(unsigned int num_groups, double sigma_t)
{
  Reset();

  num_groups_ = num_groups;
  sigma_t_.resize(num_groups, sigma_t);
  sigma_a_.resize(num_groups, sigma_t);

  ComputeDiffusionParameters();
}


//######################################################################
/**
 * Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering.
 */
void chi_physics::MultiGroupXS::
MakeSimple1(unsigned int num_groups, double sigma_t, double c)
{
  Reset();

  num_groups_ = num_groups;
  sigma_t_.resize(num_groups, sigma_t);
  transfer_matrices_.emplace_back(num_groups, num_groups);

  // When multi-group, assign half the scattering cross section
  // to within-group scattering. The other half will be used for
  // up/down-scattering.

  auto& S = transfer_matrices_.back();
  double scale = (num_groups_ == 1) ? 1.0 : 0.5;
  S.SetDiagonal(std::vector<double>(num_groups, sigma_t * c * scale));

  // Set the up/down-scattering cross sections.
  // Summary:
  //     1) The half of groups with higher energies down-scatter to the next
  //        lowest energy group half the time and to the same group half the
  //        time.
  //     2) The half of groups with lower energies less the last group
  //        down-scatter to the next lowest energy group three quarters of the
  //        time and up-scatter to the next highest energy group one quarter
  //        of the time.
  //     3) The lowest energy group has the same form as 1).

  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    //downscattering
    if (g > 0)
      S.Insert(g, g - 1, sigma_t * c * 0.5);

    //upscattering
    if (g > num_groups_ / 2)
    {
      if (g < num_groups_ - 1)
      {
        S.Insert(g, g - 1, sigma_t * c * 0.25);
        S.Insert(g, g + 1, sigma_t * c * 0.25);
      }
      else
        S.Insert(g, g - 1, sigma_t * c * 0.5);
    }
  }//for g

  ComputeAbsorption();
  ComputeDiffusionParameters();
}


//######################################################################
/** Populates the cross section from a combination of others. */
void chi_physics::MultiGroupXS::
MakeCombined(std::vector<std::pair<int, double> > &combinations)
{
  Reset();

  //pickup all xs and make sure valid
  std::vector<std::shared_ptr<chi_physics::MultiGroupXS>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; //total density of fissile materials

  //loop over cross sections
  for (auto combo : combinations)
  {
    //get the cross section from the lua stack
    std::shared_ptr<chi_physics::MultiGroupXS> xs;
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, combo.first,
                              std::string(__FUNCTION__));
    xsecs.push_back(xs);

    //increment densities
    if (xs->is_fissionable_)
    {
      is_fissionable_ = true;
      Nf_total += combo.second;
    }

    //define and check number of groups
    if (xsecs.size() == 1)
      n_grps = xs->num_groups_;
    else if (xs->num_groups_ != n_grps)
      throw std::logic_error(
          "Incompatible cross sections encountered.\n"
          "All cross sections being combined must have the "
          "same number of energy groups.");

    //increment number of precursors
    n_precs += xs->num_precursors_;
  }//for cross section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  if (Nf_total < 1.0e-12)
    is_fissionable_ = false;

  // Check to ensure that all fissionable cross sections contain either
  // prompt/delayed fission data or total fission data
  if (n_precs > 0)
    for (const auto& xs : xsecs)
      if (xs->is_fissionable_ && xs->num_precursors_ == 0)
        throw std::logic_error(
            "Incompatible cross sections encountered.\n"
            "If any fissionable cross sections specify delayed neutron "
            "data, all fissionable cross sections must specify delayed "
            "neutron data.");

  //============================================================
  // Initialize the data
  //============================================================

  num_groups_ = n_grps;
  num_precursors_ = n_precs;
  scattering_order_ = 0;
  for (const auto& xs : xsecs)
    scattering_order_ = std::max(scattering_order_,
                                 xs->scattering_order_);

  //mandatory cross sections
  sigma_t_.assign(n_grps, 0.0);
  sigma_a_.assign(n_grps, 0.0);

  //init transfer matrices only if at least one exists
  using XSPtr = chi_physics::TransportCrossSectionsPtr;
  if (std::any_of(xsecs.begin(), xsecs.end(),
                  [](const XSPtr& x)
                  { return !x->transfer_matrices_.empty(); }))
    transfer_matrices_.assign(scattering_order_ + 1,
                              chi_math::SparseMatrix(num_groups_, num_groups_));

  //init fission data
  if (is_fissionable_)
  {
    sigma_f_.assign(n_grps, 0.0);
    nu_sigma_f_.assign(n_grps, 0.0);
    production_matrix_.assign(
        num_groups_, std::vector<double>(num_groups_, 0.0));

    //init prompt/delayed fission data
    if (n_precs > 0)
    {
      nu_prompt_sigma_f_.assign(n_grps, 0.0);
      nu_delayed_sigma_f_.assign(n_grps, 0.0);
      precursors_.resize(n_precs);
    }
  }

  //============================================================
  // Combine the data
  //============================================================

  unsigned int precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    //atom density
    double N_i = combinations[x].second;

    //fraction of fissile density
    double ff_i = 0.0;
    if (xsecs[x]->is_fissionable_)
      ff_i = N_i / Nf_total;

    //============================================================
    // Combine cross sections
    //============================================================

    // Here, raw cross sections are scaled by densities and spectra by
    // fractional densities. The latter is done to preserve a unit spectra.
    for (unsigned int g = 0; g < n_grps; ++g)
    {
      sigma_t_[g] += xsecs[x]->sigma_t_[g] * N_i;
      sigma_a_[g] += xsecs[x]->sigma_a_[g] * N_i;

      if (xsecs[x]->is_fissionable_)
      {
        sigma_f_[g] += xsecs[x]->sigma_f_[g] * N_i;
        nu_sigma_f_[g] += xsecs[x]->sigma_f_[g] * N_i;
        for (unsigned int gp = 0; gp < num_groups_; ++gp)
          production_matrix_[g][gp] +=
              xsecs[x]->production_matrix_[g][gp] * N_i;

        if (n_precs > 0)
        {
          nu_prompt_sigma_f_[g] += xsecs[g]->nu_prompt_sigma_f_[g] * N_i;
          nu_delayed_sigma_f_[g] += xsecs[g]->nu_delayed_sigma_f_[g] * N_i;
        }
      }
    }//for g

    //============================================================
    // Combine precursor data
    //============================================================

    // Here, all precursors across all materials are stored. The decay
    // constants and delayed spectrum are what they are, however, some
    // special treatment must be given to the yields. Because the yield
    // tells us what fraction of delayed neutrons are produced from a
    // given family, the sum over all families must yield unity. To
    // achieve this end, we must scale all precursor yields based on
    // the fraction of the total density of materials with precursors
    // they make up.

    if (xsecs[x]->num_precursors_ > 0)
    {
      for (unsigned int j = 0; j < xsecs[x]->num_precursors_; ++j)
      {
        unsigned int count = precursor_count + j;
        const auto& precursor = xsecs[x]->precursors_[j];
        precursors_[count].decay_constant = precursor.decay_constant;
        precursors_[count].fractional_yield = precursor.fractional_yield * ff_i;
        precursors_[count].emission_spectrum = precursor.emission_spectrum;
      }//for j

      precursor_count += xsecs[x]->num_precursors_;
    }

    //============================================================
    // Set inverse velocity data
    //============================================================

    if (x == 0 && !xsecs[x]->inv_velocity_.empty())
      inv_velocity_ = xsecs[x]->inv_velocity_;
    else if (xsecs[x]->inv_velocity_ != inv_velocity_)
      throw std::logic_error(
          "Invalid cross sections encountered.\n"
          "All cross sections being combined must share a group "
          "structure. This implies that the inverse speeds for "
          "each of the cross sections must be equivalent.");

    //============================================================
    // Combine transfer matrices
    //============================================================

    // This step is somewhat tricky. The cross sections aren't guaranteed
    // to have the same sparsity patterns and therefore simply adding them
    // together has to take the sparse matrix's protection mechanisms into
    // account.

    if (!xsecs[x]->transfer_matrices_.empty())
    {
      for (unsigned int m = 0; m < xsecs[x]->scattering_order_ + 1; ++m)
      {
        auto& Sm = transfer_matrices_[m];
        const auto& Sm_other = xsecs[x]->transfer_matrices_[m];
        for (unsigned int g = 0; g < num_groups_; ++g)
        {
          const auto& cols = Sm_other.rowI_indices_[g];
          const auto& vals = Sm_other.rowI_values_[g];
          for (size_t t = 0; t < cols.size(); ++t)
            Sm.InsertAdd(g, t, vals[t] * N_i);
        }
      }
    }
  }//for cross sections

  //============================================================
  // Compute auxiliary data
  //============================================================

  ComputeDiffusionParameters();
}


//######################################################################
void chi_physics::MultiGroupXS::ComputeAbsorption()
{
  sigma_a_.assign(num_groups_, 0.0);

  // compute for a pure absorber
  if (transfer_matrices_.empty())
    for (size_t g = 0; g < num_groups_; ++g)
      sigma_a_[g] = sigma_t_[g];

  // estimate from a transfer matrix
  else
  {
    chi::log.Log0Warning()
        << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices_[0];
    for (size_t g = 0; g < num_groups_; ++g)
    {
      // estimate the scattering cross section
      double sig_s = 0.0;
      for (size_t row = 0; row < S0.NumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices_[row];
        const auto& vals = S0.rowI_values_[row];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sig_s += vals[t];
            break;
          }
      }

      sigma_a_[g] = sigma_t_[g] - sig_s;

      // TODO: Should negative absorption be allowed?
      if (sigma_a_[g] < 0.0)
        chi::log.Log0Warning()
            << "Negative absorption cross section encountered "
            << "in group " << g << " when estimating from the "
            << "transfer matrices";
    }//for g
  }//if scattering present
}


//######################################################################
/** Scale the fission data by a constant. */
void chi_physics::MultiGroupXS::ScaleFissionData(const double k)
{
  if (is_fission_scaled_)
  {
    chi::log.Log0Warning()
        << "An attempt was made to scale fission data after "
           "it had already been scaled... Nothing will be done.";
    return;
  }
  
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    nu_sigma_f_[g] /= k;
    nu_prompt_sigma_f_[g] /= k;
    nu_delayed_sigma_f_[g] /= k;

    auto& prod = production_matrix_[g];
    for (auto& val : prod)
      val /= k;
  }
  is_fission_scaled_ = true;
}


//######################################################################
chi_physics::AdjointMultiGroupXS::
AdjointMultiGroupXS(const MultiGroupXS& xs) : xs_(xs)
{
  // transpose transfer matrices
  for (unsigned int ell = 0; ell < xs_.ScatteringOrder() + 1; ++ell)
  {
    const auto& S_ell = xs_.TransferMatrix(ell);
    chi_math::SparseMatrix S_ell_transpose(xs_.NumGroups(), xs_.NumGroups());
    for (size_t g = 0; g < xs_.NumGroups(); ++g)
    {
      const size_t row_len = S_ell.rowI_indices_[g].size();
      const size_t* col_ptr = S_ell.rowI_indices_[g].data();
      const double* val_ptr = S_ell.rowI_values_[g].data();

      for (size_t j = 0; j < row_len; ++j)
        S_ell_transpose.Insert(*col_ptr++, g, *val_ptr++);

      transposed_transfer_matrices_.push_back(S_ell_transpose);
    }
  }//for ell

  // transpose production matrices
  if (xs_.IsFissionable())
  {
    const auto& F = xs_.ProductionMatrix();
    for (size_t g = 0; g < xs_.NumGroups(); ++g)
    {
      std::vector<double> F_g_transpose;
      for (size_t gp = 0; gp < xs_.NumGroups(); ++gp)
        F_g_transpose.emplace_back(F[gp][g]);
      transposed_production_matrices_.push_back(F_g_transpose);
    }
  }
}
