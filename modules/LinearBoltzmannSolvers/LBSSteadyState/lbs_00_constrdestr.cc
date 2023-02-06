#include "lbs_linear_boltzmann_solver.h"



//###################################################################
/**Constructor for LBS*/
lbs::SteadyStateSolver::SteadyStateSolver(const std::string& in_text_name) :
  chi_physics::Solver(in_text_name)
{}

/**Returns the time at which the last restart was written.*/
double lbs::SteadyStateSolver::LastRestartWrite() const
{return last_restart_write_;}

/**Returns a reference to the time at which the last restart was written.*/
double& lbs::SteadyStateSolver::LastRestartWrite()
{return last_restart_write_;}

/**Returns a reference to the solver options.*/
lbs::Options& lbs::SteadyStateSolver::Options()
{return options_;}

/**Returns the number of moments for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::SteadyStateSolver::NumMoments() const
{return num_moments_;}

/**Returns the number of groups_ for the solver. This will only be non-zero
 * after initialization.*/
size_t lbs::SteadyStateSolver::NumGroups() const
{return num_groups_;}

/**Adds a group to the list of groups_. If group id < 0, the id will be logically
 * derived from the list size. If >= 0 the id will be set to the id specified.*/
void lbs::SteadyStateSolver::AddGroup(int id)
{
  if (id < 0)
    groups_.emplace_back(static_cast<int>(groups_.size()));
  else
    groups_.emplace_back(id);
}

/**Const accessor.*/
const std::vector<lbs::LBSGroup>& lbs::SteadyStateSolver::Groups() const
{
  return groups_;
}

/**Adds a groupset to the list of groupsets_. The groupset id will be logically
 * derived from the list size.*/
void lbs::SteadyStateSolver::AddGroupset()
{
  groupsets_.emplace_back(static_cast<int>(groupsets_.size()));
}

/**Non-Const accessor.*/
std::vector<lbs::LBSGroupset>& lbs::SteadyStateSolver::Groupsets()
{
  return groupsets_;
}

/**Const accessor.*/
const std::vector<lbs::LBSGroupset>& lbs::SteadyStateSolver::Groupsets() const
{
  return groupsets_;
}

/**Obtains a reference to the spatial discretization.*/
const chi_math::SpatialDiscretization& lbs::SteadyStateSolver::
  SpatialDiscretization() const
{
  return *discretization_;
}

/**Adds a point source to the solver's point source list.*/
void lbs::SteadyStateSolver::AddPointSource(PointSource psrc)
{
  point_sources_.push_back(std::move(psrc));
}

/**Clears all the point sources from the solver's point source list.*/
void lbs::SteadyStateSolver::ClearPointSources()
{
  point_sources_.clear();
}

/**Const accessor to the list of point sources.*/
const std::vector<lbs::PointSource>& lbs::SteadyStateSolver::
  PointSources() const
{
  return point_sources_;
}

/**Obtains a reference to the unknown manager for flux-moments.*/
const chi_math::UnknownManager& lbs::SteadyStateSolver::UnknownManager() const
{
  return flux_moments_uk_man_;
}

/**Returns the local node count for the flux-moments data structures.*/
size_t lbs::SteadyStateSolver::LocalNodeCount() const
{return local_node_count_;}

/**Returns the global node count for the flux-moments data structures.*/
size_t lbs::SteadyStateSolver::GlobalNodeCount() const
{return glob_node_count_;}

/**Read/write access to source moments vector.*/
std::vector<double>& lbs::SteadyStateSolver::QMomentsLocal()
{
  return q_moments_local_;
}

/**Read/write access to exterior src moments vector.*/
std::vector<double>& lbs::SteadyStateSolver::ExtSrcMomentsLocal()
{
  return ext_src_moments_local_;
}

/**Read/write access to last updated flux vector.*/
std::vector<double>& lbs::SteadyStateSolver::PhiOldLocal()
{
  return phi_old_local_;
}

/**Read/write access to newest updated flux vector.*/
std::vector<double>& lbs::SteadyStateSolver::PhiNewLocal()
{
  return phi_new_local_;
}

/**Read/Write access to the boundary preferences.*/
std::map<uint64_t, lbs::BoundaryPreference>& lbs::SteadyStateSolver::
  BoundaryPreferences()
{
  return boundary_preferences_;
}