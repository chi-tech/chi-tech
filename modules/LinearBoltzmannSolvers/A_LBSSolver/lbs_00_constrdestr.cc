#include "lbs_solver.h"

#include "chi_log.h"

#include "ChiObject/object_maker.h"

namespace lbs
{
RegisterChiObject(lbs, LBSSolver);
}

/**Base class constructor.*/
lbs::LBSSolver::LBSSolver(const std::string& text_name)
  : chi_physics::Solver(text_name)
{
}

/**Returns the input parameters for this object.*/
chi_objects::InputParameters lbs::LBSSolver::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "LBSDatablock");

  params.AddRequiredParameter<size_t>(
    "num_groups", "The total number of groups within the solver");

  params.AddRequiredParameterArray(
    "groupsets",
    "An array of blocks each specifying the input parameters for a groupsets");

  return params;
}

/**Input parameters based construction.*/
lbs::LBSSolver::LBSSolver(const chi_objects::InputParameters& params)
  : chi_physics::Solver(params)
{
  //=================================== Make groups
  const size_t num_groups = params.GetParamValue<size_t>("num_groups");
  for (size_t g = 0; g < num_groups; ++g)
    groups_.push_back(LBSGroup(static_cast<int>(g)));

  //=================================== Make groupsets
  const auto& groupsets_array = params.GetParam("groupsets");

  const size_t num_gs = groupsets_array.NumParameters();
  for (size_t gs = 0; gs < num_gs; ++gs)
  {
    const auto& groupset_params = groupsets_array.GetParam(gs);

    chi_objects::InputParameters gs_input_params =
      LBSGroupset::GetInputParameters();
    gs_input_params.SetObjectType("LBSSolver:LBSGroupset");
    gs_input_params.AssignParameters(groupset_params);

    groupsets_.emplace_back(gs_input_params, gs, *this);
  } // for gs
}

/**Returns the source event tag used for logging the time it
 * takes to set source moments.*/
size_t lbs::LBSSolver::GetSourceEventTag() const { return source_event_tag_; }

/**Returns the time at which the last restart was written.*/
double lbs::LBSSolver::LastRestartWrite() const { return last_restart_write_; }

/**Returns a reference to the time at which the last restart was written.*/
double& lbs::LBSSolver::LastRestartWrite() { return last_restart_write_; }

/**Returns a reference to the solver options.*/
lbs::Options& lbs::LBSSolver::Options() { return options_; }

/**Returns a constant reference to the solver options.*/
const lbs::Options& lbs::LBSSolver::Options() const { return options_; }

/**Returns the number of moments for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumMoments() const { return num_moments_; }

/**Returns the number of groups for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumGroups() const { return num_groups_; }

/**Returns the number of precursors for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumPrecursors() const { return num_precursors_; }

/**Returns the maximum number of precursors, for a material, as encountered
 * accross all the materials. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::GetMaxPrecursorsPerMaterial() const
{
  return max_precursors_per_material_;
}

/**Adds a group to the list of groups. If group id < 0, the id will be logically
 * derived from the list size. If >= 0 the id will be set to the id specified.*/
void lbs::LBSSolver::AddGroup(int id)
{
  if (id < 0) groups_.emplace_back(static_cast<int>(groups_.size()));
  else
    groups_.emplace_back(id);
}

/**Const accessor.*/
const std::vector<lbs::LBSGroup>& lbs::LBSSolver::Groups() const
{
  return groups_;
}

/**Adds a groupset to the list of groupsets. The groupset id will be logically
 * derived from the list size.*/
void lbs::LBSSolver::AddGroupset()
{
  groupsets_.emplace_back(static_cast<int>(groupsets_.size()));
}

/**Non-Const accessor.*/
std::vector<lbs::LBSGroupset>& lbs::LBSSolver::Groupsets()
{
  return groupsets_;
}

/**Const accessor.*/
const std::vector<lbs::LBSGroupset>& lbs::LBSSolver::Groupsets() const
{
  return groupsets_;
}

/**Adds a point source to the solver's point source list.*/
void lbs::LBSSolver::AddPointSource(PointSource psrc)
{
  point_sources_.push_back(std::move(psrc));
}

/**Clears all the point sources from the solver's point source list.*/
void lbs::LBSSolver::ClearPointSources() { point_sources_.clear(); }

/**Const accessor to the list of point sources.*/
const std::vector<lbs::PointSource>& lbs::LBSSolver::PointSources() const
{
  return point_sources_;
}

/**Returns a reference to the map of material ids to XSs.*/
const std::map<int, lbs::XSPtr>& lbs::LBSSolver::GetMatID2XSMap() const
{
  return matid_to_xs_map_;
}

/**Returns a reference to the map of material ids to Isotropic Srcs.*/
const std::map<int, lbs::IsotropicSrcPtr>&
lbs::LBSSolver::GetMatID2IsoSrcMap() const
{
  return matid_to_src_map_;
}

/**Obtains a reference to the spatial discretization.*/
const chi_math::SpatialDiscretization&
lbs::LBSSolver::SpatialDiscretization() const
{
  return *discretization_;
}

/**Returns read-only access to the unit cell matrices.*/
const std::vector<lbs::UnitCellMatrices>&
lbs::LBSSolver::GetUnitCellMatrices() const
{
  return unit_cell_matrices_;
}

/**Obtains a reference to the grid.*/
const chi_mesh::MeshContinuum& lbs::LBSSolver::Grid() const
{
  return *grid_ptr_;
}

/**Returns a reference to the list of local cell transport views.*/
const std::vector<lbs::CellLBSView>&
lbs::LBSSolver::GetCellTransportViews() const
{
  return cell_transport_views_;
}

/**Obtains a reference to the unknown manager for flux-moments.*/
const chi_math::UnknownManager& lbs::LBSSolver::UnknownManager() const
{
  return flux_moments_uk_man_;
}

/**Returns the local node count for the flux-moments data structures.*/
size_t lbs::LBSSolver::LocalNodeCount() const { return local_node_count_; }

/**Returns the global node count for the flux-moments data structures.*/
size_t lbs::LBSSolver::GlobalNodeCount() const { return glob_node_count_; }

/**Read/write access to source moments vector.*/
std::vector<double>& lbs::LBSSolver::QMomentsLocal()
{
  return q_moments_local_;
}
/**Read access to source moments vector.*/
const std::vector<double>& lbs::LBSSolver::QMomentsLocal() const
{
  return q_moments_local_;
}

/**Read/write access to exterior src moments vector.*/
std::vector<double>& lbs::LBSSolver::ExtSrcMomentsLocal()
{
  return ext_src_moments_local_;
}

/**Read access to exterior src moments vector.*/
const std::vector<double>& lbs::LBSSolver::ExtSrcMomentsLocal() const
{
  return ext_src_moments_local_;
}

/**Read/write access to last updated flux vector.*/
std::vector<double>& lbs::LBSSolver::PhiOldLocal() { return phi_old_local_; }

/**Read access to last updated flux vector.*/
const std::vector<double>& lbs::LBSSolver::PhiOldLocal() const
{
  return phi_old_local_;
}

/**Read/write access to newest updated flux vector.*/
std::vector<double>& lbs::LBSSolver::PhiNewLocal() { return phi_new_local_; }
/**Read access to newest updated flux vector.*/
const std::vector<double>& lbs::LBSSolver::PhiNewLocal() const
{
  return phi_new_local_;
}

/**Read/write access to newest updated precursors vector.*/
std::vector<double>& lbs::LBSSolver::PrecursorsNewLocal()
{
  return phi_new_local_;
}
/**Read access to newest updated precursors vector.*/
const std::vector<double>& lbs::LBSSolver::PrecursorsNewLocal() const
{
  return phi_new_local_;
}

/**Read/write access to newest updated angular flux vector.*/
std::vector<VecDbl>& lbs::LBSSolver::PsiNewLocal() { return psi_new_local_; }
/**Read access to newest updated angular flux vector.*/
const std::vector<VecDbl>& lbs::LBSSolver::PsiNewLocal() const
{
  return psi_new_local_;
}

/**Returns the sweep boundaries as a read only reference*/
const std::map<uint64_t, std::shared_ptr<SweepBndry>>&
lbs::LBSSolver::SweepBoundaries() const
{
  return sweep_boundaries_;
}

/**Read/Write access to the boundary preferences.*/
std::map<uint64_t, lbs::BoundaryPreference>&
lbs::LBSSolver::BoundaryPreferences()
{
  return boundary_preferences_;
}

/**Gets the local and global number of iterative unknowns. This normally is
 * only the flux moments, however, the sweep based solvers might include
 * delayed angular fluxes in this number.*/
std::pair<size_t, size_t> lbs::LBSSolver::GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  return {num_local_phi_dofs, num_globl_phi_dofs};
}

/**Gets the local handle of a flux-moment based field function.*/
size_t lbs::LBSSolver::MapPhiFieldFunction(size_t g, size_t m) const
{
  ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
                    std::string("Failure to map phi field function g") +
                      std::to_string(g) + " m" + std::to_string(m));

  return phi_field_functions_local_map_.at({g, m});
}

/**Returns the local handle to the power generation field function, if
 * enabled.*/
size_t lbs::LBSSolver::GetHandleToPowerGenFieldFunc() const
{
  ChiLogicalErrorIf(not options_.power_field_function_on,
                    "Called when options_.power_field_function_on == false");

  return power_gen_fieldfunc_local_handle_;
}