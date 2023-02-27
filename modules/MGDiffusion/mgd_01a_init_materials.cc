#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include <algorithm>
//============================================= assemble matrix A
void mg_diffusion::Solver::Initialize_Materials(std::set<int>& material_ids)
{
  chi::log.Log0Verbose1() << "Initializing Materials";

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  const size_t num_physics_mats = chi::material_stack.size();
  bool first_material_read = true;

  for (const int& mat_id : material_ids)
  {
    auto current_material = chi::GetStackItemPtr(chi::material_stack,
                                                 mat_id, __FUNCTION__);
    materials_list << "Material id " << mat_id;

    //====================================== Check valid ids
    if (mat_id < 0) {throw std::logic_error(
        "MG-diff-InitMaterials: Cells encountered with no assigned material.");}
    if (static_cast<size_t>(mat_id) >= num_physics_mats) {
      throw std::logic_error("MG-diff-InitMaterials: Cells encountered with "
                             "material id that matches no material in physics material library.");}

    //====================================== Extract properties
    using MatProperty = chi_physics::PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties_)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          std::static_pointer_cast<chi_physics::TransportCrossSections>(property);
        matid_to_xs_map[mat_id] = transp_xs;
        found_transport_xs = true;
        if (first_material_read)
          mg_diffusion::Solver::num_groups_ = transp_xs->num_groups_;

      }//transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source =
          std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g_.size() < mg_diffusion::Solver::num_groups_)
        {
          chi::log.LogAllWarning()
            << "MG-Diff-InitMaterials: Isotropic Multigroup source specified in "
            << "material \"" << current_material->name_ << "\" has fewer "
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
        << "MG-Diff-InitMaterials: Found no transport cross-section property for "
        << "material \"" << current_material->name_ << "\".";
      chi::Exit(EXIT_FAILURE);
    }
    //====================================== Check number of groups legal
    if (matid_to_xs_map[mat_id]->num_groups_ != mg_diffusion::Solver::num_groups_)
    {
      chi::log.LogAllError()
          << "MG-Diff-InitMaterials: Found material \"" << current_material->name_ << "\" has "
          << matid_to_xs_map[mat_id]->num_groups_ << " groups and"
          << " the simulation has " << mg_diffusion::Solver::num_groups_ << " groups."
          << " The material must have the same number of groups.";
      chi::Exit(EXIT_FAILURE);
    }

    //====================================== Check number of moments
    if (matid_to_xs_map[mat_id]->scattering_order_ > 1)
    {
      chi::log.Log0Warning()
          << "MG-Diff-InitMaterials: Found material \"" << current_material->name_ << "\" has "
          << "a scattering order of "
          << matid_to_xs_map[mat_id]->scattering_order_ << " and"
          << " the simulation has a scattering order of One (MG-Diff)"
          << " The higher moments will therefore not be used.";
    }

    materials_list
        << " number of moments "
        << matid_to_xs_map[mat_id]->transfer_matrices_.size() << "\n";

    first_material_read = false;
  }//for material id


  chi::log.Log()
    << "Materials Initialized:\n" << materials_list.str() << "\n";

  MPI_Barrier(MPI_COMM_WORLD);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize Diffusion
  //                                                   properties
  chi::log.Log() << "Computing diffusion parameters.";

  for (const auto& mat_id_xs : matid_to_xs_map)
    mat_id_xs.second->ComputeDiffusionParameters();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute last fast group
  // initialize last fast group
  chi::log.Log() << "Computing last fast group.";
  unsigned int lfg = mg_diffusion::Solver::num_groups_;

  if (mg_diffusion::Solver::num_groups_ > 1)
  {
    // loop over all materials
    for (const auto &mat_id_xs: matid_to_xs_map)
    {
      // get the P0 transfer matrix
      const auto &S = mat_id_xs.second->transfer_matrices_[0];
      // loop over all row of the transfer matrix
      const int G = static_cast<int>(mg_diffusion::Solver::num_groups_);
      for (int g = G-1; g >=0 ; --g)
      {
        for (const auto &[row_g, gp, sigma_sm]: S.Row(g))
        {
          if ((std::fabs(sigma_sm) > 1e-10) && (gp > row_g))
            lfg = std::min(lfg, static_cast<unsigned int>(row_g));
        }
      }
    }// loop on materials
  }//if num_groups>1

  mg_diffusion::Solver::last_fast_group_ = lfg;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute two-grid params
  do_two_grid_ = basic_options_("do_two_grid").BoolValue();
  if ((lfg == num_groups_) && do_two_grid_)
  {
    chi::log.Log0Error() << "Two-grid is not possible with no upscattering.";
    do_two_grid_ = false;
    chi::Exit(EXIT_FAILURE);
  }
  if (do_two_grid_)
  {
    chi::log.Log() << "Compute_TwoGrid_Params";
    mg_diffusion::Solver::Compute_TwoGrid_Params();
  }

}
