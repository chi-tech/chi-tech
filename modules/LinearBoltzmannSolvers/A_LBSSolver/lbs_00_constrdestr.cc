#include "lbs_solver.h"

/**Base class constructor.*/
lbs::LBSSolver::LBSSolver(const std::string& text_name) :
  chi_physics::Solver(text_name)
{

}

/**Returns the source event tag used for logging the time it
 * takes to set source moments.*/
size_t lbs::LBSSolver::GetSourceEventTag() const
{
  return source_event_tag_;
}

/**Returns the time at which the last restart was written.*/
double lbs::LBSSolver::LastRestartWrite() const
{return last_restart_write_;}

/**Returns a reference to the time at which the last restart was written.*/
double& lbs::LBSSolver::LastRestartWrite()
{return last_restart_write_;}

/**Returns a reference to the solver options.*/
lbs::Options& lbs::LBSSolver::Options()
{return options_;}

/**Returns a constant reference to the solver options.*/
const lbs::Options& lbs::LBSSolver::Options() const
{return options_;}

/**Returns the number of moments for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumMoments() const
{return num_moments_;}

/**Returns the number of groups for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumGroups() const
{return num_groups_;}

/**Returns the number of precursors for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::NumPrecursors() const
{return num_precursors_;}

/**Returns the maximum number of precursors, for a material, as encountered
 * accross all the materials. This will only be non-zero
 * after initialization.*/
size_t lbs::LBSSolver::GetMaxPrecursorsPerMaterial() const
{return max_precursors_per_material_;}

/**Adds a group to the list of groups. If group id < 0, the id will be logically
 * derived from the list size. If >= 0 the id will be set to the id specified.*/
void lbs::LBSSolver::AddGroup(int id)
{
  if (id < 0)
    groups_.emplace_back(static_cast<int>(groups_.size()));
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
void lbs::LBSSolver::ClearPointSources()
{
  point_sources_.clear();
}

/**Const accessor to the list of point sources.*/
const std::vector<lbs::PointSource>& lbs::LBSSolver::
PointSources() const
{
  return point_sources_;
}

/**Returns a reference to the map of material ids to XSs.*/
const std::map<int, lbs::XSPtr>& lbs::LBSSolver::GetMatID2XSMap() const
{
  return matid_to_xs_map_;
}

/**Returns a reference to the map of material ids to Isotropic Srcs.*/
const std::map<int, lbs::IsotropicSrcPtr>& lbs::LBSSolver::GetMatID2IsoSrcMap() const
{
  return matid_to_src_map_;
}

/**Obtains a reference to the spatial discretization.*/
const chi_math::SpatialDiscretization& lbs::LBSSolver::
SpatialDiscretization() const
{
  return *discretization_;
}

/**Obtains a reference to the grid.*/
const chi_mesh::MeshContinuum& lbs::LBSSolver::
  Grid() const
{
  return *grid_ptr_;
}

/**Returns a reference to the list of local cell transport views.*/
const std::vector<lbs::CellLBSView>& lbs::LBSSolver::GetCellTransportViews() const
{
  return cell_transport_views_;
}

/**Obtains a reference to the unknown manager for flux-moments.*/
const chi_math::UnknownManager& lbs::LBSSolver::UnknownManager() const
{
  return flux_moments_uk_man_;
}

/**Returns the local node count for the flux-moments data structures.*/
size_t lbs::LBSSolver::LocalNodeCount() const
{return local_node_count_;}

/**Returns the global node count for the flux-moments data structures.*/
size_t lbs::LBSSolver::GlobalNodeCount() const
{return glob_node_count_;}

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
std::vector<double>& lbs::LBSSolver::PhiOldLocal()
{
  return phi_old_local_;
}

/**Read access to last updated flux vector.*/
const std::vector<double>& lbs::LBSSolver::PhiOldLocal() const
{
  return phi_old_local_;
}

/**Read/write access to newest updated flux vector.*/
std::vector<double>& lbs::LBSSolver::PhiNewLocal()
{
  return phi_new_local_;
}
/**Read access to newest updated flux vector.*/
const std::vector<double>& lbs::LBSSolver::PhiNewLocal() const
{
  return phi_new_local_;
}

/**Read/Write access to the boundary preferences.*/
std::map<uint64_t, lbs::BoundaryPreference>& lbs::LBSSolver::
BoundaryPreferences()
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