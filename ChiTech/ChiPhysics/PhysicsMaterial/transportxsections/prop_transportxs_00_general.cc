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

  sigma_tg.clear();
  sigma_tg.resize(in_G,in_sigmat);
  sigma_fg.resize(in_G,0.0);
  sigma_captg.resize(in_G,0.0);
  chi_g.resize(in_G,0.0);
  nu_sigma_fg.resize(in_G,0.0);
  nu_p_sigma_fg.resize(in_G,0.0);
  nu_d_sigma_fg.resize(in_G,0.0);
  ddt_coeff.resize(in_G,0.0);

  transfer_matrix.emplace_back(in_G,in_G);
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

  sigma_tg.resize(in_G,in_sigmat);
  sigma_tg.clear();
  sigma_tg.resize(in_G,in_sigmat);
  sigma_fg.resize(in_G,0.0);
  sigma_captg.resize(in_G,0.0);
  chi_g.resize(in_G,0.0);
  nu_sigma_fg.resize(in_G,0.0);
  nu_p_sigma_fg.resize(in_G,0.0);
  nu_d_sigma_fg.resize(in_G,0.0);
  ddt_coeff.resize(in_G,0.0);

  transfer_matrix.emplace_back(in_G,in_G);

  auto& ref_matrix = transfer_matrix.back();

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
  int num_grps_G = 0;
  int num_precursors_J = 0;
  int count = 0;
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
    ++count;
  }//for auto
  if (Nf_total < 1.0e-28)
    Nf_total = 1.0; //Avoids divide by 0 when non-fissile

  //======================================== Combine 1D cross-sections
  this->num_groups = num_grps_G;
  this->num_precursors = num_precursors_J;
  sigma_tg.clear();
  sigma_fg.clear();
  sigma_captg.clear();
  chi_g.clear();
  nu_sigma_fg.clear();
  nu_p_sigma_fg.clear();
  nu_d_sigma_fg.clear();
  ddt_coeff.clear();
  lambda.clear();
  gamma.clear();
  chi_d.clear();

  sigma_tg.resize(num_grps_G,0.0);
  sigma_fg.resize(num_grps_G,0.0);
  sigma_captg.resize(num_grps_G,0.0);
  chi_g.resize(num_grps_G,0.0);
  nu_sigma_fg.resize(num_grps_G,0.0);
  nu_p_sigma_fg.resize(num_grps_G,0.0);
  nu_d_sigma_fg.resize(num_grps_G,0.0);
  ddt_coeff.resize(num_grps_G,0.0);
  lambda.resize(num_precursors_J,0.0);
  gamma.resize(num_precursors_J,0.0);
  chi_d.resize(num_grps_G);
  for (int g=0; g < num_groups; ++g)
    chi_d[g].resize(num_precursors_J,0.0);
  precursor_map.resize(num_precursors_J,0);

  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    this->scattering_order = std::max(this->scattering_order, cross_secs[x]->scattering_order);

    double N_i = combinations[x].second;
    double f_i = N_i/N_total;
    double ff_i = N_i/Nf_total;

    for (int g=0; g<num_grps_G; g++)
    {
      sigma_tg     [g] += cross_secs[x]->sigma_tg     [g] * N_i;
      sigma_fg     [g] += cross_secs[x]->sigma_fg     [g] * N_i;
      sigma_captg  [g] += cross_secs[x]->sigma_captg  [g] * N_i;
      chi_g        [g] += cross_secs[x]->chi_g        [g] * ff_i;
      nu_sigma_fg  [g] += cross_secs[x]->nu_sigma_fg  [g] * N_i;
      nu_p_sigma_fg[g] += cross_secs[x]->nu_p_sigma_fg[g] * N_i;
      nu_d_sigma_fg[g] += cross_secs[x]->nu_d_sigma_fg[g] * N_i;
      ddt_coeff    [g] += cross_secs[x]->ddt_coeff    [g] * f_i;
    }
    if ((cross_secs[x]->is_fissile) and (cross_secs[x]->num_precursors > 0))
    {
      for (int j=0; j<num_precursors_J; ++j)
      {
        lambda[j] += cross_secs[x]->lambda[j] * ff_i;
        gamma [j] += cross_secs[x]->gamma [j] * ff_i;
        for (int g=0; g < num_groups; g++)
          chi_d[g][j] += cross_secs[x]->chi_d[g][j] * ff_i;
      }
    }
  }

  //======================================== Combine transfer matrices
  // This step is somewhat tricky. The cross-sections
  // aren't guaranteed to have the same sparsity patterns
  // and therefore simply adding them together has to take
  // the sparse matrix's protection mechanisms into account.
  transfer_matrix.clear();
  transfer_matrix.resize(this->scattering_order + 1,
                         chi_math::SparseMatrix(num_grps_G,num_grps_G));
  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    for (int m=0; m<(cross_secs[x]->scattering_order + 1); ++m)
    {
      auto& xs_tm = cross_secs[x]->transfer_matrix[m];
      for (int i=0; i < num_groups; ++i)
      {
        for (auto j : xs_tm.rowI_indices[i])
        {
          double value = xs_tm.ValueIJ(i,j)*combinations[x].second;
          transfer_matrix[m].InsertAdd(i,j,value);
        }
      }//for i
    }//for m
  }//for xs
}