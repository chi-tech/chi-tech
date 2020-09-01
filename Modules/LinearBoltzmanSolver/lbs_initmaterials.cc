#include "lbs_linear_boltzman_solver.h"

#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>


extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <algorithm>

//###################################################################
/**Initializes default materials and physics materials.*/
void LinearBoltzman::Solver::InitMaterials(std::set<int>& material_ids)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "Initializing Materials";

  //================================================== Determine dsa flags
  bool develop_wgdsa = false;
  for (int gs=0; gs<group_sets.size(); gs++)
  {
    if ((group_sets[gs]->apply_wgdsa) ||
        (group_sets[gs]->apply_tgdsa))
      develop_wgdsa = true;
  }

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check materials found
  int num_physics_mats = chi_physics_handler.material_stack.size();
  matid_to_xs_map.resize(num_physics_mats,-1);
  matid_to_src_map.resize(num_physics_mats,-1);
  std::set<int>::iterator matID;
  for (matID  = material_ids.begin();
       matID != material_ids.end();
       matID++)
  {
    int id = *matID;

    materials_list << "Material id " << *matID;

    //====================================== Check valid ids
    if (id<0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Cells encountered with no assigned material.";
      exit(EXIT_FAILURE);
    }
    if (id>=num_physics_mats)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Cells encountered with material id that matches no"
           "material in physics material library.";
      exit(EXIT_FAILURE);
    }

    //====================================== Extract properties
    using MatProperty = chi_physics::PropertyType;
    chi_physics::Material* cur_mat = chi_physics_handler.material_stack[id];
    int num_props = cur_mat->properties.size();
    bool found_transport_xs = false;
    for (int p=0; p<num_props; p++)
    {
      if (cur_mat->properties[p]->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          (chi_physics::TransportCrossSections*)cur_mat->properties[p];
        material_xs.push_back(transp_xs);
        matid_to_xs_map[id] = material_xs.size()-1;
        found_transport_xs = true;
      }
      if (cur_mat->properties[p]->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source =
          (chi_physics::IsotropicMultiGrpSource*)cur_mat->properties[p];

        if (mg_source->source_value_g.size() < groups.size())
        {
          chi_log.Log(LOG_ALLWARNING)
            << "LBS-InitMaterials: Isotropic Multigroup source specified in "
            << "material \"" << cur_mat->name << "\" has fewer energy groups than "
            << "called for in the simulation. Source will be ignored.";
          //exit(EXIT_FAILURE);
        }
        else
        {
          material_srcs.push_back(mg_source);
          matid_to_src_map[id] = material_srcs.size()-1;
        }
      }
    }//for property

    //====================================== Check valid property
    if (!found_transport_xs)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Found no transport cross-section property for "
        << "material \"" << cur_mat->name << "\".";
      exit(EXIT_FAILURE);
    }

    //====================================== Check number of groups legal
    if (material_xs[matid_to_xs_map[id]]->G < groups.size())
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Found material \"" << cur_mat->name << "\" has "
        << material_xs[matid_to_xs_map[id]]->G << " groups and"
        << " the simulation has " << groups.size() << " groups."
        << " The material must have a greater or equal amount of groups.";
      exit(EXIT_FAILURE);
    }

    //====================================== Check number of moments
    if (material_xs[matid_to_xs_map[id]]->L < options.scattering_order)
    {
      chi_log.Log(LOG_0WARNING)
        << "LBS-InitMaterials: Found material \"" << cur_mat->name << "\" has "
        << "a scattering order of "
        << material_xs[matid_to_xs_map[id]]->L << " and"
        << " the simulation has a scattering order of "
        << options.scattering_order << "."
        << " The higher moments will therefore not be used.";
    }

    materials_list
    << " number of moments "
    << material_xs[matid_to_xs_map[id]]->transfer_matrix.size() << "\n";
  }//for material id

  chi_log.Log(LOG_0)
  << "Materials Initialized:\n"
  << materials_list.str() << "\n";

  MPI_Barrier(MPI_COMM_WORLD);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize WGDSA stuff
  if (develop_wgdsa)
  {
    chi_log.Log(LOG_0) << "Computing diffusion parameters.";

    for (int xs=0; xs<material_xs.size(); xs++)
    {
      material_xs[xs]->ComputeDiffusionParameters();
      auto rel_location = std::find(matid_to_xs_map.begin(),
                                    matid_to_xs_map.end(), xs);
      size_t mat_id = rel_location - matid_to_xs_map.begin();
      if (rel_location != matid_to_xs_map.end())
        chi_log.Log(LOG_0VERBOSE_1) << "Material " << mat_id;
    }
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do more stuff

}
