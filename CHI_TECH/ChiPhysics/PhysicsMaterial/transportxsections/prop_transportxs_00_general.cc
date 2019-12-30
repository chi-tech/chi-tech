#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics chi_physics_handler;

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**Default constructor.*/
chi_physics::TransportCrossSections::TransportCrossSections() :
  chi_physics::MaterialProperty(PropertyType::TRANSPORT_XSECTIONS)
{
  G = 0;
  L = 0;

  diffusion_initialized = false;
  scattering_initialized = false;
}

//###################################################################
/**Makes a simple material with no transfer matrix just sigma_t.*/
void chi_physics::TransportCrossSections::
  MakeSimple0(int in_G, double in_sigmat)
{
  G = in_G;
  sigma_tg.clear();
  sigma_tg.resize(in_G,in_sigmat);
  sigma_fg.resize(in_G,0.0);
  sigma_captg.resize(in_G,0.0);
  chi_g.resize(in_G,0.0);
  nu_sigma_fg.resize(in_G,0.0);

  transfer_matrix.push_back(chi_math::SparseMatrix(in_G,in_G));
}

//###################################################################
/**Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering.*/
void chi_physics::TransportCrossSections::
  MakeSimple1(int in_G, double in_sigmat, double c)
{
  G = in_G;

  sigma_tg.resize(in_G,in_sigmat);
  sigma_tg.clear();
  sigma_tg.resize(in_G,in_sigmat);
  sigma_fg.resize(in_G,0.0);
  sigma_captg.resize(in_G,0.0);
  chi_g.resize(in_G,0.0);
  nu_sigma_fg.resize(in_G,0.0);

  transfer_matrix.push_back(chi_math::SparseMatrix(in_G,in_G));

  if (G == 1)
    transfer_matrix[0].SetDiagonal(std::vector<double>(in_G,in_sigmat*c));
  else
    transfer_matrix[0].SetDiagonal(std::vector<double>(in_G,in_sigmat*c*2.0/4.0));

  for (int g=0; g<in_G; g++)
  {
    //Downscattering
    if (g>0)
    {
      transfer_matrix[0][g][g-1] = in_sigmat*c*2.0/4.0;
    }


    //Upscattering
    if (g>(in_G/2))
    {
      if (g<(in_G-1))
      {
        transfer_matrix[0][g][g-1] = in_sigmat*c*1.0/4.0;
        transfer_matrix[0][g][g+1] = in_sigmat*c*1.0/4.0;
      }
      else
      {
        transfer_matrix[0][g][g-1] = in_sigmat*c*2.0/4.0;
      }
    }

  }

//  chi_log.Log(LOG_0WARNING) << transfer_matrix[0].PrintS();
}

//###################################################################
/**Populates the cross-section from a combination of others.*/
void chi_physics::TransportCrossSections::
  MakeCombined(std::vector<std::pair<int, double> > &combinations)
{
  //======================================== Pickup all xs and make sure valid
  std::vector<chi_physics::TransportCrossSections*> cross_secs;
  cross_secs.reserve(combinations.size());
  int num_grps_G=0;
  int count=0;
  for (auto& combo : combinations)
  {
    chi_physics::TransportCrossSections* xs;
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

    //============================ Check number of groups
    if (cross_secs.size() == 1)
      num_grps_G = xs->G;
    else
      if (cross_secs[count-1]->G != num_grps_G)
        chi_log.Log(LOG_ALLERROR)
          << "In call to TransportCrossSections::MakeCombined: "
          << "all cross-sections must have the same number of groups.";
    ++count;
  }

  //======================================== Combine 1D cross-sections
  this->G = num_grps_G;
  sigma_tg.clear();
  sigma_fg.clear();
  sigma_captg.clear();
  chi_g.clear();
  nu_sigma_fg.clear();

  sigma_tg.resize(num_grps_G,0.0);
  sigma_fg.resize(num_grps_G,0.0);
  sigma_captg.resize(num_grps_G,0.0);
  chi_g.resize(num_grps_G,0.0);
  nu_sigma_fg.resize(num_grps_G,0.0);
  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    this->L = std::max(this->L,cross_secs[x]->L);
    for (int g=0; g<G; g++)
    {
      sigma_tg   [g] += cross_secs[x]->sigma_tg   [g]*combinations[x].second;
      sigma_fg   [g] += cross_secs[x]->sigma_fg   [g]*combinations[x].second;
      sigma_captg[g] += cross_secs[x]->sigma_captg[g]*combinations[x].second;
      chi_g      [g] += cross_secs[x]->chi_g      [g]*combinations[x].second;
      nu_sigma_fg[g] += cross_secs[x]->nu_sigma_fg[g]*combinations[x].second;
    }
  }

  //======================================== Combine transfer matrices
  // This step is somewhat tricky. The cross-sections
  // aren't guaranteed to have the same sparsity patterns
  // and therefore simply adding them together has to take
  // the sparse matrix's protection mechanisms into account.
  transfer_matrix.clear();
  transfer_matrix.resize(this->L+1,
                         chi_math::SparseMatrix(num_grps_G,num_grps_G));
  for (size_t x=0; x<cross_secs.size(); ++x)
  {
    for (int m=0; m<(cross_secs[x]->L+1); ++m)
    {
      auto& xs_tm = cross_secs[x]->transfer_matrix[m];
      for (int i=0; i<G; ++i)
      {
        for (auto j : xs_tm.inds_rowI[i])
        {
          double value = xs_tm.ValueIJ(i,j)*combinations[x].second;
          transfer_matrix[m].InsertAdd(i,j,value);
        }
      }//for i
    }//for m
  }//for xs

}