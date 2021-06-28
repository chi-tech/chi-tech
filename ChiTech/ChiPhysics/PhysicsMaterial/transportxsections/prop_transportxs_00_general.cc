#include "material_property_transportxsections.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

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
  double N_total = 0.0;
  double Nf_total = 0.0;
  for (auto combo : combinations)
  {
    std::shared_ptr<chi_physics::TransportCrossSections> xs;
    try {
      xs = chi_physics_handler.trnsprt_xs_stack.at(combo.first);
    }
    catch(const std::out_of_range& o){
      chi_log.Log(LOG_ALLERROR)
        << "ERROR: Invalid cross-section handle"
        << " in call to chiPhysicsMaterialSetProperty."
        << std::endl;
      exit(EXIT_FAILURE);
    }
    
    cross_secs.push_back(xs);

    // Increment combo factor totals
    N_total += combo.second;
    if (xs->is_fissile) {
      this->is_fissile = true;
      Nf_total += combo.second;
    }

    //============================ Check number of groups
    if (cross_secs.size() == 1)
    {
      num_grps_G = xs->num_groups;
      if (xs->is_fissile)
        num_precursors_J = xs->num_precursors;
    }
    else
    {
      if (xs->num_groups != num_grps_G)
      {
        chi_log.Log(LOG_ALLERROR)
          << "In call to TransportCrossSections::MakeCombined: "
          << "all cross-sections must have the same number of groups.";
        exit(EXIT_FAILURE);
      }
      if (xs->is_fissile)
      {
        if (num_precursors_J == 0) 
          num_precursors_J = xs->num_precursors;
        else
        {
          if (xs->num_precursors != num_precursors_J)
          {
            chi_log.Log(LOG_ALLERROR)
              << "In call to TransportCrossSections::MakeCombined: "
              << "all fissile cross-sections must have the same number "
              << "of precursors.";
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }//for auto
  if (Nf_total < 1.0e-28)
    Nf_total = 1.0; //Avoids divide by 0 when non-fissile

  //======================================== Combine 1D cross-sections
  this->num_groups = num_grps_G;
  this->num_precursors = num_precursors_J;
  sigma_t.clear();
  sigma_f.clear();
  sigma_a.clear();
  chi.clear();
  chi_prompt.clear();
  nu.clear();
  nu_prompt.clear();
  nu_delayed.clear();
  nu_sigma_f.clear();
  nu_prompt_sigma_f.clear();
  nu_delayed_sigma_f.clear();
  inv_velocity.clear();
  precursor_lambda.clear();
  precursor_yield.clear();
  chi_delayed.clear();

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
  for (int g=0; g < num_groups; ++g)
    chi_delayed[g].resize(num_precursors_J, 0.0);

  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    scattering_order = std::max(this->scattering_order,
                                cross_secs[x]->scattering_order);

    double N_i = combinations[x].second;
    double f_i = N_i/N_total;
    double ff_i = N_i/Nf_total;

    for (int g=0; g<num_grps_G; g++)
    {
      sigma_t     [g] += cross_secs[x]->sigma_t     [g] * N_i;
      sigma_f     [g] += cross_secs[x]->sigma_f     [g] * N_i;
      sigma_a     [g] += cross_secs[x]->sigma_a     [g] * N_i;
      chi        [g] += cross_secs[x]->chi        [g] * ff_i;
      chi_prompt [g] += cross_secs[x]->chi_prompt[g] * ff_i;
      nu           [g] += cross_secs[x]->nu           [g] * ff_i;
      nu_prompt    [g] += cross_secs[x]->nu_prompt    [g] * ff_i;
      nu_delayed   [g] += cross_secs[x]->nu_delayed   [g] * ff_i;
      nu_sigma_f  [g] += cross_secs[x]->nu_sigma_f  [g] * N_i;
      nu_prompt_sigma_f[g] += cross_secs[x]->nu_prompt_sigma_f[g] * N_i;
      nu_delayed_sigma_f[g] += cross_secs[x]->nu_delayed_sigma_f[g] * N_i;
      inv_velocity    [g] += cross_secs[x]->inv_velocity    [g] * f_i;
    }
    if ((cross_secs[x]->is_fissile) and (cross_secs[x]->num_precursors > 0))
    {
      for (int j=0; j<num_precursors_J; ++j)
      {
        precursor_lambda[j] += cross_secs[x]->precursor_lambda[j] * ff_i;
        precursor_yield [j] += cross_secs[x]->precursor_yield [j] * ff_i;
        for (int g=0; g < num_groups; g++)
          chi_delayed[g][j] += cross_secs[x]->chi_delayed[g][j] * ff_i;
      }
    }
  }

  //======================================== Combine transfer matrices
  // This step is somewhat tricky. The cross-sections
  // aren't guaranteed to have the same sparsity patterns
  // and therefore simply adding them together has to take
  // the sparse matrix's protection mechanisms into account.
  transfer_matrices.clear();
  transfer_matrices.resize(this->scattering_order + 1,
                           chi_math::SparseMatrix(num_grps_G,num_grps_G));
  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    for (int m=0; m<(cross_secs[x]->scattering_order + 1); ++m)
    {
      auto& xs_tm = cross_secs[x]->transfer_matrices[m];
      for (int i=0; i < num_groups; ++i)
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