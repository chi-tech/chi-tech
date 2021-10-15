#include "lbs_linear_boltzmann_solver.h"

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
void LinearBoltzmann::Solver::InitMaterials(std::set<int>& material_ids)
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing Materials";

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  int num_physics_mats = chi_physics_handler.material_stack.size();
  material_xs.clear();
  material_srcs.clear();
  matid_to_xs_map.assign(num_physics_mats,-1);
  matid_to_src_map.assign(num_physics_mats,-1);
  for (const int& mat_id : material_ids)
  {
    auto current_material = chi_physics_handler.material_stack[mat_id];
    materials_list << "Material id " << mat_id;

    //====================================== Check valid ids
    if (mat_id < 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Cells encountered with no assigned material.";
      exit(EXIT_FAILURE);
    }
    if (mat_id >= num_physics_mats)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Cells encountered with material id "
        << "that matches no material in physics material library.";
      exit(EXIT_FAILURE);
    }

    //====================================== Extract properties
    using MatProperty = chi_physics::PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          std::static_pointer_cast<chi_physics::TransportCrossSections>(property);
        material_xs.push_back(transp_xs);
        matid_to_xs_map[mat_id] = material_xs.size() - 1;
        found_transport_xs = true;
      }//transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source =
          std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g.size() < groups.size())
        {
          chi_log.Log(LOG_ALLWARNING)
            << "LBS-InitMaterials: Isotropic Multigroup source specified in "
            << "material \"" << current_material->name << "\" has fewer "
            << "energy groups than called for in the simulation. "
            << "Source will be ignored.";
        }
        else
        {
          material_srcs.push_back(mg_source);
          matid_to_src_map[mat_id] = material_srcs.size() - 1;
        }
      }//P0 source
    }//for property

    //====================================== Check valid property
    if (!found_transport_xs)
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Found no transport cross-section property for "
        << "material \"" << current_material->name << "\".";
      exit(EXIT_FAILURE);
    }

    //====================================== Check number of groups legal
    if (material_xs[matid_to_xs_map[mat_id]]->num_groups < groups.size())
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBS-InitMaterials: Found material \"" << current_material->name << "\" has "
        << material_xs[matid_to_xs_map[mat_id]]->num_groups << " groups and"
        << " the simulation has " << groups.size() << " groups."
        << " The material must have a greater or equal amount of groups.";
      exit(EXIT_FAILURE);
    }

    //====================================== Check number of moments
    if (material_xs[matid_to_xs_map[mat_id]]->scattering_order < options.scattering_order)
    {
      chi_log.Log(LOG_0WARNING)
        << "LBS-InitMaterials: Found material \"" << current_material->name << "\" has "
        << "a scattering order of "
        << material_xs[matid_to_xs_map[mat_id]]->scattering_order << " and"
        << " the simulation has a scattering order of "
        << options.scattering_order << "."
        << " The higher moments will therefore not be used.";
    }
    
    materials_list
      << " number of moments "
      << material_xs[matid_to_xs_map[mat_id]]->transfer_matrices.size() << "\n";
  }//for material id

  num_groups = groups.size();

  chi_log.Log(LOG_0)
    << "Materials Initialized:\n" << materials_list.str() << "\n";

  MPI_Barrier(MPI_COMM_WORLD);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize precursor
  //                                                   properties
  num_precursors = 0;
  max_precursors_per_material = 0;
  for (auto& xs : material_xs)
  {
    num_precursors += xs->num_precursors;
    if (xs->num_precursors > max_precursors_per_material)
      max_precursors_per_material = xs->num_precursors;
  }
  if (num_precursors == 0)
    options.use_precursors = false;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize Diffusion
  //                                                   properties
  bool develop_diffusion_properties = false;
  for (auto& group_set : groupsets)
  {
    if (group_set.apply_wgdsa || group_set.apply_tgdsa)
      develop_diffusion_properties = true;
  }

  if (develop_diffusion_properties)
  {
    chi_log.Log(LOG_0) << "Computing diffusion parameters.";

    for (auto& xs : material_xs)
      xs->ComputeDiffusionParameters();
  }

}
