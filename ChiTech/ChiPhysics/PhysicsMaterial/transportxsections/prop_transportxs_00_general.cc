#include "chi_runtime.h"
#include "material_property_transportxsections.h"

#include "chi_runtime.h"
#include "chi_log.h"
;

//###################################################################
/**Default constructor.*/
chi_physics::TransportCrossSections::TransportCrossSections() :
  chi_physics::MaterialProperty(PropertyType::TRANSPORT_XSECTIONS)
{
  num_groups = 0;
  scattering_order = 0;

  diffusion_initialized = false;
  scattering_initialized = false;
}

//###################################################################
/**Makes a simple material with no transfer matrix just sigma_t.*/
void chi_physics::TransportCrossSections::
  MakeSimple0(int in_G, double in_sigmat)
{
  //======================================== Clear any previous data
  Reset();

  num_groups = in_G;
  scattering_order = 0;

  sigma_t.clear();
  sigma_t.resize(in_G, in_sigmat);
  sigma_f.resize(in_G, 0.0);
  sigma_a.resize(in_G, 0.0);
  chi.resize(in_G, 0.0);
  chi_prompt.resize(in_G, 0.0);
  nu_sigma_f.resize(in_G, 0.0);
  nu_prompt_sigma_f.resize(in_G, 0.0);
  nu_delayed_sigma_f.resize(in_G, 0.0);
  inv_velocity.resize(in_G, 0.0);

  transfer_matrices.emplace_back(in_G, in_G);
}

//###################################################################
/**Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering.*/
void chi_physics::TransportCrossSections::
  MakeSimple1(int in_G, double in_sigmat, double c)
{
  //======================================== Clear any previous data
  Reset();

  num_groups = in_G;
  scattering_order = 0;

  sigma_t.resize(in_G, in_sigmat);
  sigma_f.resize(in_G, 0.0);
  sigma_a.resize(in_G, 0.0);
  chi.resize(in_G, 0.0);
  chi_prompt.resize(in_G, 0.0);
  nu_sigma_f.resize(in_G, 0.0);
  nu_prompt_sigma_f.resize(in_G, 0.0);
  nu_delayed_sigma_f.resize(in_G, 0.0);
  inv_velocity.resize(in_G, 0.0);

  transfer_matrices.emplace_back(in_G, in_G);

  auto& ref_matrix = transfer_matrices.back();

  if (num_groups == 1)
    ref_matrix.SetDiagonal(std::vector<double>(in_G,in_sigmat*c));
  else
    ref_matrix.SetDiagonal(std::vector<double>(in_G,in_sigmat*c*2.0/4.0));

  for (int g=0; g<in_G; g++)
  {
    //Downscattering
    if (g>0)
      ref_matrix.Insert(g,g-1,in_sigmat*c*2.0/4.0);

    //Upscattering
    if (g>(in_G/2))
    {
      if (g<(in_G-1))
      {
        ref_matrix.Insert(g,g-1,in_sigmat*c*1.0/4.0);
        ref_matrix.Insert(g,g+1,in_sigmat*c*1.0/4.0);
      }
      else
        ref_matrix.Insert(g,g-1,in_sigmat*c*2.0/4.0);
    }
  }//for g

  sigma_a = ComputeAbsorptionXSFromTransfer();
}

//###################################################################
/**Populates the cross-section from a combination of others.*/
void chi_physics::TransportCrossSections::
  MakeCombined(std::vector<std::pair<int, double> > &combinations)
{
  //======================================== Clear any previous data
  Reset();

  //======================================== Pickup all xs and make sure valid
  std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> cross_secs;
  cross_secs.reserve(combinations.size());
  size_t num_grps_G = 0;
  size_t num_precursors_J = 0;

  double N_total  = 0.0; // total density
  double Nf_total = 0.0; // total density of fissile materials
  double Np_total = 0.0; // total density of materials with precursors

  for (auto combo : combinations)
  {
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
    
    cross_secs.push_back(xs);

    // Increment combo factor totals
    N_total += combo.second;
    if (xs->is_fissile)
    {
      this->is_fissile = true;
      Nf_total += combo.second;

      if (xs->num_precursors > 0)
        Np_total += combo.second;
    }

    //============================ Define and check number of groups
    if (cross_secs.size() == 1)
      num_grps_G = xs->num_groups;
    else if (xs->num_groups != num_grps_G)
    {
      chi::log.LogAllError()
        << "In call to " << __FUNCTION__ << ": "
        << "all cross-sections must have the same number of groups.";
      chi::Exit(EXIT_FAILURE);
    }

    //============================ Increment number of precursors
    if (not xs->is_fissile and xs->num_precursors > 0)
    {
      chi::log.LogAllError()
          << "In call to " << __FUNCTION__ << ": "
          << "only fissile materials are allowed to have delayed "
          << "neutron precursors.";
      chi::Exit(EXIT_FAILURE);
    }
    num_precursors_J += xs->num_precursors;
  }//for cross-section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  double eps = std::numeric_limits<double>::epsilon(); //machine precision
  if (Nf_total < eps)
    this->is_fissile = false;
  if (Np_total < eps)
    this->num_precursors = 0;

  //======================================== Combine 1D cross-sections
  this->num_groups = num_grps_G;
  this->num_precursors = num_precursors_J;

  sigma_t.resize(num_grps_G, 0.0);
  sigma_f.resize(num_grps_G, 0.0);
  sigma_a.resize(num_grps_G, 0.0);
  chi.resize(num_grps_G, 0.0);
  chi_prompt.resize(num_grps_G, 0.0);
  nu.resize(num_grps_G,0.0);
  nu_prompt.resize(num_grps_G,0.0);
  nu_delayed.resize(num_grps_G,0.0);
  nu_sigma_f.resize(num_grps_G, 0.0);
  nu_prompt_sigma_f.resize(num_grps_G, 0.0);
  nu_delayed_sigma_f.resize(num_grps_G, 0.0);
  inv_velocity.resize(num_grps_G, 0.0);
  precursor_lambda.resize(num_precursors_J, 0.0);
  precursor_yield.resize(num_precursors_J, 0.0);
  chi_delayed.resize(num_grps_G);
  for (size_t g = 0; g < num_groups; ++g)
    chi_delayed[g].resize(num_precursors_J, 0.0);

  size_t precursor_count = 0;
  for (size_t x = 0; x < cross_secs.size(); ++x)
  {
    scattering_order = std::max(this->scattering_order,
                                cross_secs[x]->scattering_order);

    // Atom density
    double N_i = combinations[x].second;

    // Fraction of fissile density
    double ff_i = 0.0;
    if (cross_secs[x]->is_fissile)
      ff_i = N_i / Nf_total;

    // Fraction of precursor density
    double pf_i = 0.0;
    if (cross_secs[x]->num_precursors > 0)
      pf_i = N_i / Np_total;

    //======================================== Combine cross-sections
    for (size_t g = 0; g < num_grps_G; ++g)
    {
      sigma_t     [g] += cross_secs[x]->sigma_t     [g] * N_i;
      sigma_f     [g] += cross_secs[x]->sigma_f     [g] * N_i;
      sigma_a     [g] += cross_secs[x]->sigma_a     [g] * N_i;

      chi         [g] += cross_secs[x]->chi        [g] * ff_i;
      chi_prompt  [g] += cross_secs[x]->chi_prompt [g] * ff_i;

      nu          [g] += cross_secs[x]->nu           [g] * ff_i;
      nu_prompt   [g] += cross_secs[x]->nu_prompt    [g] * ff_i;
      nu_delayed  [g] += cross_secs[x]->nu_delayed   [g] * ff_i;

      nu_sigma_f        [g] += cross_secs[x]->nu_sigma_f        [g] * N_i;
      nu_prompt_sigma_f [g] += cross_secs[x]->nu_prompt_sigma_f [g] * N_i;
      nu_delayed_sigma_f[g] += cross_secs[x]->nu_delayed_sigma_f[g] * N_i;

      if (x == 0)
        inv_velocity[g] = cross_secs[x]->inv_velocity[g];
      else if (inv_velocity[g] != cross_secs[x]->inv_velocity[g])
        chi::log.LogAllWarning()
            << "In call to " << __FUNCTION__
            << ": all materials must have the same inv_velocity "
            << "term per group. Invalid value encountered in "
            << "material " << x << " group " << g << ". Using the "
            << "value from the first cross-section set.";
    }

    //======================================== Compute precursor
    // Here, all precursors across all materials are stored.
    // The decay constants and delayed spectrum are what they are,
    // however, some special treatment must be given to the yields.
    // Because the yield tells us what fraction of delayed neutrons
    // are produced from a given family, the sum over all families
    // must yield unity. To achieve this end, we must scale all
    // precursor yields based on the fraction of the total density
    // of materials with precursors they make up.
    if (cross_secs[x]->num_precursors > 0)
    {
      for (size_t j = 0; j < cross_secs[x]->num_precursors; ++j)
      {
        size_t count = precursor_count + j;
        precursor_lambda[count] = cross_secs[x]->precursor_lambda[j];
        precursor_yield [count] = cross_secs[x]->precursor_yield [j] * pf_i;
        for (size_t g = 0; g < num_groups; ++g)
          chi_delayed[g][count] = cross_secs[x]->chi_delayed[g][j];
      }
      precursor_count += cross_secs[x]->num_precursors;
    }
  }//for cross sections

  //======================================== Combine transfer matrices
  // This step is somewhat tricky. The cross-sections
  // aren't guaranteed to have the same sparsity patterns
  // and therefore simply adding them together has to take
  // the sparse matrix's protection mechanisms into account.
  transfer_matrices.clear();
  transfer_matrices.resize(this->scattering_order + 1,
                           chi_math::SparseMatrix(num_grps_G,num_grps_G));
  for (size_t x = 0; x < cross_secs.size(); ++x)
  {
    for (int m=0; m<(cross_secs[x]->scattering_order + 1); ++m)
    {
      auto& xs_tm = cross_secs[x]->transfer_matrices[m];
      for (size_t i = 0; i < num_groups; ++i)
      {
        for (auto j : xs_tm.rowI_indices[i])
        {
          double value = xs_tm.ValueIJ(i,j)*combinations[x].second;
          transfer_matrices[m].InsertAdd(i, j, value);
        }
      }//for i
    }//for m
  }//for xs
}

std::vector<chi_physics::TransportCrossSections::GrpVal>
  chi_physics::TransportCrossSections::ComputeAbsorptionXSFromTransfer()
{
  auto GetMatrixColumnSum = [](const chi_math::SparseMatrix& matrix, size_t col)
  {
    double sum = 0.0;

    for (size_t row=0; row<matrix.NumRows(); ++row)
    {
      const auto& row_col_indices = matrix.rowI_indices[row];
      const auto& row_values = matrix.rowI_values[row];
      for (size_t j=0; j<row_col_indices.size(); ++j)
        if (row_col_indices[j] == col)
          sum += row_values[j];
    }

    return sum;
  };

  std::vector<GrpVal> temp_sigma_a(num_groups, 0.0);
  if (not transfer_matrices.empty())
    for (size_t g = 0; g < num_groups; ++g)
      temp_sigma_a[g] = sigma_t[g] - sigma_f[g] -
        GetMatrixColumnSum(transfer_matrices[0], g);

  return temp_sigma_a;
}