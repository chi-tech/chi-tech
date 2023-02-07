#include "gs_context.h"

#include <petscksp.h>

#include "LinearBoltzmannSolvers/LBSSteadyState/lbs_linear_boltzmann_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

#define sc_double static_cast<double>

namespace lbs
{

template<>
std::pair<int64_t, int64_t> GSContext<Mat,Vec>::SystemSize()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();
  const size_t num_moments      = lbs_solver_.NumMoments();

  const size_t groupset_numgrps = groupset_.groups.size();
  const auto num_delayed_psi_info = groupset_.angle_agg.GetNumDelayedAngularDOFs();
  const size_t local_size = local_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.first;
  const size_t globl_size = globl_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.second;
  const size_t num_angles = groupset_.quadrature->abscissae.size();
  const size_t num_psi_global = globl_node_count *
                                num_angles *
                                groupset_.groups.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info_)
  {
    chi::log.Log()
      << "Total number of angular unknowns: "
      << num_psi_global
      << "\n"
      << "Number of lagged angular unknowns: "
      << num_delayed_psi_globl << "("
      << std::setprecision(2)
      << sc_double(num_delayed_psi_globl)*100 / sc_double(num_psi_global)
      << "%)";
  }

  return {static_cast<int64_t>(local_size),
          static_cast<int64_t>(globl_size)};
}

}