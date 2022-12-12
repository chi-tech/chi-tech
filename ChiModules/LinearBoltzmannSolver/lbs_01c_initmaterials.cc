#include "lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"


#include <algorithm>

//###################################################################
/**Initializes default materials and physics materials.*/
void lbs::SteadySolver::InitMaterials(std::set<int>& material_ids)
{
  chi::log.Log0Verbose1() << "Initializing Materials";

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  const size_t num_physics_mats = chi::material_stack.size();

  for (const int& mat_id : material_ids)
  {
    auto current_material = chi::GetStackItemPtr(chi::material_stack,
                                                 mat_id, __FUNCTION__);
    materials_list << "Material id " << mat_id;

    //====================================== Check valid ids
    if (mat_id < 0) {throw std::logic_error(
      "LBS-InitMaterials: Cells encountered with no assigned material.");}
    if (static_cast<size_t>(mat_id) >= num_physics_mats) {
      throw std::logic_error("LBS-InitMaterials: Cells encountered with "
        "material id that matches no material in physics material library.");}

    //====================================== Extract properties
    using MatProperty = chi_physics::PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          std::static_pointer_cast<chi_physics::TransportCrossSections>(property);
        matid_to_xs_map[mat_id] = transp_xs;
        found_transport_xs = true;
      }//transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source =
          std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g.size() < groups.size())
        {
          chi::log.LogAllWarning()
            << "LBS-InitMaterials: Isotropic Multigroup source specified in "
            << "material \"" << current_material->name << "\" has fewer "
            << "energy groups than called for in the simulation. "
            << "Source will be ignored.";
        }
        else
        {
          matid_to_src_map[mat_id] = mg_source;
        }
      }//P0 source
    }//for property

    //====================================== Check valid property
    if (!found_transport_xs)
    {
      chi::log.LogAllError()
        << "LBS-InitMaterials: Found no transport cross-section property for "
        << "material \"" << current_material->name << "\".";
      chi::Exit(EXIT_FAILURE);
    }
    //====================================== Check number of groups legal
    if (matid_to_xs_map[mat_id]->num_groups < groups.size())
    {
      chi::log.LogAllError()
        << "LBS-InitMaterials: Found material \"" << current_material->name << "\" has "
        << matid_to_xs_map[mat_id]->num_groups << " groups and"
        << " the simulation has " << groups.size() << " groups."
        << " The material must have a greater or equal amount of groups.";
      chi::Exit(EXIT_FAILURE);
    }

    //====================================== Check number of moments
    if (matid_to_xs_map[mat_id]->scattering_order < options.scattering_order)
    {
      chi::log.Log0Warning()
        << "LBS-InitMaterials: Found material \"" << current_material->name << "\" has "
        << "a scattering order of "
        << matid_to_xs_map[mat_id]->scattering_order << " and"
        << " the simulation has a scattering order of "
        << options.scattering_order << "."
        << " The higher moments will therefore not be used.";
    }

    materials_list
      << " number of moments "
      << matid_to_xs_map[mat_id]->transfer_matrices.size() << "\n";
  }//for material id

  num_groups = groups.size();

  chi::log.Log()
    << "Materials Initialized:\n" << materials_list.str() << "\n";

  MPI_Barrier(MPI_COMM_WORLD);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize precursor
  //                                                   properties
  num_precursors = 0;
  max_precursors_per_material = 0;
  for (const auto& mat_id_xs : matid_to_xs_map)
  {
    const auto& xs = mat_id_xs.second;
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
    chi::log.Log() << "Computing diffusion parameters.";

    for (const auto& mat_id_xs : matid_to_xs_map)
      mat_id_xs.second->ComputeDiffusionParameters();
  }

}
