#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>


//######################################################################
/** Default constructor. */
chi_physics::TransportCrossSections::TransportCrossSections() :
  chi_physics::MaterialProperty(PropertyType::TRANSPORT_XSECTIONS)
{
  num_groups = 0;
  scattering_order = 0;

  diffusion_initialized = false;
  scattering_initialized = false;
}


//######################################################################
void chi_physics::TransportCrossSections::
Reset()
{
  num_groups = 0;
  scattering_order = 0;
  num_precursors = 0;
  is_fissionable = false;

  sigma_t.clear();
  sigma_f.clear();
  sigma_a.clear();

  nu_sigma_f.clear();
  nu_prompt_sigma_f.clear();
  nu_delayed_sigma_f.clear();

  inv_velocity.clear();

  transfer_matrices.clear();
  production_matrix.clear();

  precursors.clear();

  //Diffusion quantities
  diffusion_initialized = false;
  diffusion_coeff.clear();
  sigma_removal.clear();
  sigma_s_gtog.clear();

  //Monte-Carlo quantities
  scattering_initialized = false;
  cdf_gprime_g.clear();
  scat_angles_gprime_g.clear();
}


//######################################################################
/** Makes a simple material with no transfer matrix just sigma_t. */
void chi_physics::TransportCrossSections::
MakeSimple0(int n_grps, double sigma)
{
  Reset();

  num_groups = n_grps;
  sigma_t.resize(n_grps, sigma);
  sigma_a.resize(n_grps, sigma);

  ComputeDiffusionParameters();
}


//######################################################################
/**
 * Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering.
 */
void chi_physics::TransportCrossSections::
MakeSimple1(int n_grps, double sigma, double c)
{
  Reset();

  num_groups = n_grps;
  sigma_t.resize(n_grps, sigma);
  transfer_matrices.emplace_back(n_grps, n_grps);

  // When multi-group, assign half the scattering cross section
  // to within-group scattering. The other half will be used for
  // up/down-scattering.

  auto& S = transfer_matrices.back();
  double scale = (num_groups == 1)? 1.0 : 0.5;
  S.SetDiagonal(std::vector<double>(n_grps, sigma * c * scale));

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

  for (unsigned int g = 0; g < num_groups; ++g)
  {
    //downscattering
    if (g > 0)
      S.Insert(g, g - 1, sigma * c * 0.5);

    //upscattering
    if (g > num_groups / 2)
    {
      if (g < num_groups - 1)
      {
        S.Insert(g, g - 1, sigma * c * 0.25);
        S.Insert(g, g + 1, sigma * c * 0.25);
      }
      else
        S.Insert(g, g - 1, sigma * c * 0.5);
    }
  }//for g

  ComputeAbsorption();
  ComputeDiffusionParameters();
}


//######################################################################
/** Populates the cross section from a combination of others. */
void chi_physics::TransportCrossSections::
MakeCombined(std::vector<std::pair<int, double> > &combinations)
{
  Reset();

  //pickup all xs and make sure valid
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; //total density of fissile materials

  //loop over cross sections
  for (auto combo : combinations)
  {
    //get the cross section from the lua stack
    std::shared_ptr<chi_physics::TransportCrossSections> xs;
    xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, combo.first,
                              std::string(__FUNCTION__));
    xsecs.push_back(xs);

    //increment densities
    if (xs->is_fissionable)
    {
      is_fissionable = true;
      Nf_total += combo.second;
    }

    //define and check number of groups
    if (xsecs.size() == 1)
      n_grps = xs->num_groups;
    else if (xs->num_groups != n_grps)
      throw std::logic_error(
          "Incompatible cross sections encountered.\n"
          "All cross sections being combined must have the "
          "same number of energy groups.");

    //increment number of precursors
    n_precs += xs->num_precursors;
  }//for cross section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  if (Nf_total < 1.0e-12)
    is_fissionable = false;

  // Check to ensure that all fissionable cross sections contain either
  // prompt/delayed fission data or total fission data
  if (n_precs > 0)
    for (const auto& xs : xsecs)
      if (xs->is_fissionable && xs->num_precursors == 0)
        throw std::logic_error(
            "Incompatible cross sections encountered.\n"
            "If any fissionable cross sections specify delayed neutron "
            "data, all fissionable cross sections must specify delayed "
            "neutron data.");

  //============================================================
  // Initialize the data
  //============================================================

  num_groups = n_grps;
  num_precursors = n_precs;
  scattering_order = 0;
  for (const auto& xs : xsecs)
    scattering_order = std::max(scattering_order,
                                xs->scattering_order);

  //mandatory cross sections
  sigma_t.assign(n_grps, 0.0);
  sigma_a.assign(n_grps, 0.0);

  //init transfer matrices only if at least one exists
  using XSPtr = chi_physics::TransportCrossSectionsPtr;
  if (std::any_of(xsecs.begin(), xsecs.end(),
                  [](const XSPtr& x)
                  { return !x->transfer_matrices.empty(); }))
    transfer_matrices.assign(scattering_order + 1,
                             TransferMatrix(num_groups, num_groups));

  //init fission data
  if (is_fissionable)
  {
    sigma_f.assign(n_grps, 0.0);
    nu_sigma_f.assign(n_grps, 0.0);
    production_matrix.assign(
        num_groups, std::vector<double>(num_groups, 0.0));

    //init prompt/delayed fission data
    if (n_precs > 0)
    {
      nu_prompt_sigma_f.assign(n_grps, 0.0);
      nu_delayed_sigma_f.assign(n_grps, 0.0);
      precursors.resize(n_precs);
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
    if (xsecs[x]->is_fissionable)
      ff_i = N_i / Nf_total;

    //============================================================
    // Combine cross sections
    //============================================================

    // Here, raw cross sections are scaled by densities and spectra by
    // fractional densities. The latter is done to preserve a unit spectra.
    for (unsigned int g = 0; g < n_grps; ++g)
    {
      sigma_t[g] += xsecs[x]->sigma_t[g] * N_i;
      sigma_a[g] += xsecs[x]->sigma_a[g] * N_i;

      if (xsecs[x]->is_fissionable)
      {
        sigma_f[g] += xsecs[x]->sigma_f[g] * N_i;
        nu_sigma_f[g] += xsecs[x]->sigma_f[g] * N_i;
        for (unsigned int gp = 0; gp < num_groups; ++gp)
          production_matrix[g][gp] +=
              xsecs[x]->production_matrix[g][gp] * N_i;

        if (n_precs > 0)
        {
          nu_prompt_sigma_f[g] += xsecs[g]->nu_prompt_sigma_f[g] * N_i;
          nu_delayed_sigma_f[g] += xsecs[g]->nu_delayed_sigma_f[g] * N_i;
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

    if (xsecs[x]->num_precursors > 0)
    {
      for (unsigned int j = 0; j < xsecs[x]->num_precursors; ++j)
      {
        unsigned int count = precursor_count + j;
        const auto& precursor = xsecs[x]->precursors[j];
        precursors[count].decay_constant = precursor.decay_constant;
        precursors[count].fractional_yield = precursor.fractional_yield * ff_i;
        precursors[count].emission_spectrum = precursor.emission_spectrum;
      }//for j

      precursor_count += xsecs[x]->num_precursors;
    }

    //============================================================
    // Set inverse velocity data
    //============================================================

    if (x == 0 && !xsecs[x]->inv_velocity.empty())
      inv_velocity = xsecs[x]->inv_velocity;
    else if (xsecs[x]->inv_velocity != inv_velocity)
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

    if (!xsecs[x]->transfer_matrices.empty())
    {
      for (unsigned int m = 0; m < xsecs[x]->scattering_order + 1; ++m)
      {
        auto& Sm = transfer_matrices[m];
        const auto& Sm_other = xsecs[x]->transfer_matrices[m];
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          const auto& cols = Sm_other.rowI_indices[g];
          const auto& vals = Sm_other.rowI_values[g];
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
void chi_physics::TransportCrossSections::
ComputeAbsorption()
{
  sigma_a.assign(num_groups, 0.0);

  // compute for a pure absorber
  if (transfer_matrices.empty())
    for (size_t g = 0; g < num_groups; ++g)
      sigma_a[g] = sigma_t[g];

  // estimate from a transfer matrix
  else
  {
    chi::log.Log0Warning()
        << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices[0];
    for (size_t g = 0; g < num_groups; ++g)
    {
      // estimate the scattering cross section
      double sig_s = 0.0;
      for (size_t row = 0; row < S0.NumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices[row];
        const auto& vals = S0.rowI_values[row];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sig_s += vals[t];
            break;
          }
      }

      sigma_a[g] = sigma_t[g] - sig_s;

      // TODO: Should negative absorption be allowed?
      if (sigma_a[g] < 0.0)
        chi::log.Log0Warning()
            << "Negative absorption cross section encountered "
            << "in group " << g << " when estimating from the "
            << "transfer matrices";
    }//for g
  }//if scattering present
}


//######################################################################
/**Scale the fission data by a constant.*/
void chi_physics::TransportCrossSections::
ScaleFissionData(const double k)
{
  if (is_fission_scaled)
  {
    chi::log.Log0Warning()
        << "An attempt was made to scale fission data after "
           "it had already been scaled... Nothing will be done.";
    return;
  }
  
  for (unsigned int g = 0; g < num_groups; ++g)
  {
    nu_sigma_f[g] /= k;
    nu_prompt_sigma_f[g] /= k;
    nu_delayed_sigma_f[g] /= k;

    auto& prod = production_matrix[g];
    for (auto& val : prod)
      val /= k;
  }
  is_fission_scaled = true;
}
