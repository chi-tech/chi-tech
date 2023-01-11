#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"


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

  chi.clear();
  chi_prompt.clear();
  chi_delayed.clear();

  nu.clear();
  nu_prompt.clear();
  nu_delayed.clear();

  nu_sigma_f.clear();
  nu_prompt_sigma_f.clear();
  nu_delayed_sigma_f.clear();

  precursor_lambda.clear();
  precursor_yield.clear();

  inv_velocity.clear();

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
  //clear any previous data
  Reset();

  //define the cross-section data
  num_groups = n_grps;

  sigma_t.resize(n_grps, sigma);
  sigma_a.resize(n_grps, sigma);
  sigma_f.resize(n_grps, 0.0);

  chi.resize(n_grps, 0.0);
  chi_prompt.resize(n_grps, 0.0);

  nu_sigma_f.resize(n_grps, 0.0);
  nu_prompt_sigma_f.resize(n_grps, 0.0);
  nu_delayed_sigma_f.resize(n_grps, 0.0);

  inv_velocity.resize(n_grps, 0.0);

  transfer_matrices.emplace_back(n_grps, n_grps);

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
  //clear any previous data
  Reset();

  //define the cross-section data
  num_groups = n_grps;
  scattering_order = 0;

  sigma_t.resize(n_grps, sigma);
  sigma_a.resize(n_grps, 0.0);
  sigma_f.resize(n_grps, 0.0);

  chi.resize(n_grps, 0.0);
  chi_prompt.resize(n_grps, 0.0);

  nu_sigma_f.resize(n_grps, 0.0);
  nu_prompt_sigma_f.resize(n_grps, 0.0);
  nu_delayed_sigma_f.resize(n_grps, 0.0);

  inv_velocity.resize(n_grps, 0.0);

  transfer_matrices.emplace_back(n_grps, n_grps);

  auto& S = transfer_matrices.back();
  double scale = (num_groups == 1)? 1.0 : 0.5;
  S.SetDiagonal(std::vector<double>(n_grps, sigma * c * scale));

  for (unsigned int g = 0; g < num_groups; ++g)
  {
    //Downscattering
    if (g > 0)
      S.Insert(g, g - 1, sigma * c * 0.5);

    //Upscattering
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
/** Populates the cross-section from a combination of others. */
void chi_physics::TransportCrossSections::
MakeCombined(std::vector<std::pair<int, double> > &combinations)
{
  //clear any previous data
  Reset();

  //pickup all xs and make sure valid
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;

  double Nf_total = 0.0; //total density of fissile materials
  double Np_total = 0.0; //total density of materials with precursors

  //loop over cross-sections
  for (auto combo : combinations)
  {
    //get the cross-section from the lua stack
    std::shared_ptr<chi_physics::TransportCrossSections> xs;
    try
    {
      xs = chi::GetStackItemPtr(chi::trnsprt_xs_stack, combo.first);
    }
    catch(const std::out_of_range& o)
    {
      chi::log.LogAllError()
        << "ERROR: Invalid cross-section handle"
        << " in call to chiPhysicsMaterialSetProperty."
        << std::endl;
      chi::Exit(EXIT_FAILURE);
    }
    xsecs.push_back(xs);

    //increment densities
    if (xs->is_fissionable)
    {
      this->is_fissionable = true;
      Nf_total += combo.second;

      if (xs->num_precursors > 0)
        Np_total += combo.second;
    }

    //define and check number of groups
    if (xsecs.size() == 1)
      n_grps = xs->num_groups;
    else if (xs->num_groups != n_grps)
    {
      chi::log.LogAllError()
        << "In call to " << __FUNCTION__ << ": "
        << "all cross-sections must have the same number of groups.";
      chi::Exit(EXIT_FAILURE);
    }

    //increment number of precursors
    if (!xs->is_fissionable && xs->num_precursors > 0)
    {
      chi::log.LogAllError()
          << "In call to " << __FUNCTION__ << ": "
          << "Only fissionable materials are allowed to have delayed "
          << "neutron precursors.";
      chi::Exit(EXIT_FAILURE);
    }
    n_precs += xs->num_precursors;
  }//for cross-section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.

  double eps = 1.0e-12;
  if (Nf_total < eps)
    is_fissionable = false;
  if (Np_total < eps)
    num_precursors = 0;

  //============================================================
  // Initialize the data
  //============================================================

  num_groups = n_grps;
  num_precursors = n_precs;

  sigma_t.resize(n_grps, 0.0);
  sigma_f.resize(n_grps, 0.0);
  sigma_a.resize(n_grps, 0.0);
  chi.resize(n_grps, 0.0);
  chi_prompt.resize(n_grps, 0.0);
  nu.resize(n_grps,0.0);
  nu_prompt.resize(n_grps,0.0);
  nu_delayed.resize(n_grps,0.0);
  beta.resize(n_grps, 0.0);
  nu_sigma_f.resize(n_grps, 0.0);
  nu_prompt_sigma_f.resize(n_grps, 0.0);
  nu_delayed_sigma_f.resize(n_grps, 0.0);
  inv_velocity.resize(n_grps, 0.0);
  precursor_lambda.resize(n_precs, 0.0);
  precursor_yield.resize(n_precs, 0.0);
  chi_delayed.resize(n_grps, std::vector<double>(n_precs, 0.0));

  unsigned int precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    scattering_order = std::max(scattering_order,
                                xsecs[x]->scattering_order);

    //atom density
    double N_i = combinations[x].second;

    //fraction of fissile density
    double ff_i = 0.0;
    if (xsecs[x]->is_fissionable)
      ff_i = N_i / Nf_total;

    //fraction of precursor density
    double pf_i = 0.0;
    if (xsecs[x]->num_precursors > 0)
      pf_i = N_i / Np_total;

    //============================================================
    // Combine cross-sections
    //============================================================

    // Here, raw cross-sections are scaled by densities and
    // spectra by fractional densities. The latter is done to
    // preserve unit spectra. The inverse velocity term must be
    // the same across all cross-section sets, so a check is
    // performed to ensure this is the case.

    for (unsigned int g = 0; g < n_grps; ++g)
    {
      sigma_t[g] += xsecs[x]->sigma_t[g] * N_i;
      sigma_a[g] += xsecs[x]->sigma_a[g] * N_i;
      sigma_f[g] += xsecs[x]->sigma_f[g] * N_i;

      chi[g] += xsecs[x]->chi[g] * ff_i;
      chi_prompt[g] += xsecs[x]->chi_prompt[g] * ff_i;

      nu[g] += xsecs[x]->nu[g] * ff_i;
      nu_prompt[g] += xsecs[x]->nu_prompt[g] * ff_i;
      nu_delayed[g] += xsecs[x]->nu_delayed[g] * ff_i;
      beta[g] += xsecs[x]->beta[g] * ff_i;

      nu_sigma_f[g] += xsecs[x]->nu_sigma_f[g] * N_i;
      nu_prompt_sigma_f[g] += xsecs[x]->nu_prompt_sigma_f[g] * N_i;
      nu_delayed_sigma_f[g] += xsecs[x]->nu_delayed_sigma_f[g] * N_i;

      if (x == 0)
        inv_velocity[g] = xsecs[x]->inv_velocity[g];
      else if (inv_velocity[g] != xsecs[x]->inv_velocity[g])
      {
        chi::log.LogAllError()
            << "In call to " << __FUNCTION__ << ": "
            << "All materials must have the same inv_velocity "
            << "term per group. Invalid value encountered in "
            << "material " << x << " group " << g << ". Using the "
            << "value from the first cross-section set.";
        chi::Exit(EXIT_FAILURE);
      }
    }

    //============================================================
    // Compute precursors
    //============================================================

    // Here, all precursors across all materials are stored.
    // The decay constants and delayed spectrum are what they are,
    // however, some special treatment must be given to the yields.
    // Because the yield tells us what fraction of delayed neutrons
    // are produced from a given family, the sum over all families
    // must yield unity. To achieve this end, we must scale all
    // precursor yields based on the fraction of the total density
    // of materials with precursors they make up.

    if (xsecs[x]->num_precursors > 0)
    {
      for (unsigned int j = 0; j < xsecs[x]->num_precursors; ++j)
      {
        unsigned int count = precursor_count + j;
        precursor_lambda[count] = xsecs[x]->precursor_lambda[j];
        precursor_yield [count] = xsecs[x]->precursor_yield [j] * pf_i;
        for (size_t g = 0; g < num_groups; ++g)
          chi_delayed[g][count] = xsecs[x]->chi_delayed[g][j];
      }
      precursor_count += xsecs[x]->num_precursors;
    }
  }//for cross sections

  //============================================================
  // Combine transfer matrices
  //============================================================

  // This step is somewhat tricky. The cross-sections
  // aren't guaranteed to have the same sparsity patterns
  // and therefore simply adding them together has to take
  // the sparse matrix's protection mechanisms into account.

  transfer_matrices.clear();
  transfer_matrices.resize(scattering_order + 1,
                           chi_math::SparseMatrix(n_grps, n_grps));
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    for (int m=0; m<(xsecs[x]->scattering_order + 1); ++m)
    {
      auto& xs_tm = xsecs[x]->transfer_matrices[m];
      for (size_t i = 0; i < num_groups; ++i)
      {
        for (auto j : xs_tm.rowI_indices[i])
        {
          double value = xs_tm.ValueIJ(i,j) * combinations[x].second;
          transfer_matrices[m].InsertAdd(i, j, value);
        }
      }//for i
    }//for m
  }//for xs

  Finalize();
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
    chi::log.LogAllWarning()
        << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices[0];
    for (size_t g = 0; g < num_groups; ++g)
    {
      // estimate the scattering cross-section
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

      // TODO: Decide whether this should be a warning or an error.
      if (sigma_a[g] < 0.0)
        chi::log.LogAllWarning()
            << "Negative absorption cross-section encountered "
            << "in group " << g << " when estimating from the "
            << "transfer matrices";
    }//for g
  }//if scattering present

}