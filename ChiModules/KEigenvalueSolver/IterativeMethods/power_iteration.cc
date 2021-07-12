#include "../k_eigenvalue_solver.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "chi_log.h"

extern ChiLog& chi_log;

#include "chi_mpi.h"

extern ChiMPI& chi_mpi;

#include "ChiTimer/chi_timer.h"

extern ChiTimer chi_program_timer;

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

using namespace LinearBoltzmann;

#include <iomanip>

//###################################################################
/**Power iterative scheme for k-eigenvalue calculations.
 * Note that this routine currently only works when the problem
 * is defined by a single groupset.
*/
void KEigenvalue::Solver::PowerIteration()
{
  chi_log.Log(LOG_0)
      << "\n\n********** Solving k-eigenvalue problem with "
      << "the Power Method.\n\n";

  LBSGroupset& groupset = group_sets[0];

  groupset.angle_agg.ZeroIncomingDelayedPsi();

  //======================================== Setup sweep chunk
  auto sweep_chunk = SetSweepChunk(groupset);
  MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                     groupset.angle_agg,
                                     *sweep_chunk);

  //======================================== Tool the sweep chunk
  sweep_scheduler.sweep_chunk.SetDestinationPhi(phi_new_local);

  //======================================== Initial guess
  phi_prev_local.assign(phi_prev_local.size(), 1.0);
  ScopedCopySTLvectors(groupset, phi_prev_local, phi_old_local);

  //======================================== Start power iterations
  double F_prev = 1.0;
  double k_eff_prev = 1.0;
  int nit = 0;      //number of iterations
  bool converged = false;
  while (nit < max_iterations)
  {

    //============================== Clear source moments
    q_moments_local.assign(q_moments_local.size(), 0.0);

    //============================== Set the fission source
    SetKSource(groupset, q_moments_local,
               APPLY_AGS_FISSION_SOURCE |
               APPLY_WGS_FISSION_SOURCE);

    //============================== Converge the scattering source with
    //                               a fixed fission source
    if (groupset.iterative_method == IterativeMethod::CLASSICRICHARDSON)
    {
      ClassicRichardson(groupset, sweep_scheduler,
                        APPLY_WGS_SCATTER_SOURCE |
                        APPLY_AGS_SCATTER_SOURCE,
                        options.verbose_inner_iterations);
    }
    else if (groupset.iterative_method == IterativeMethod::GMRES)
    {
      GMRES(groupset, sweep_scheduler,
            APPLY_WGS_SCATTER_SOURCE,           //lhs_scope
            APPLY_AGS_SCATTER_SOURCE,           //rhs_scope
            options.verbose_inner_iterations);
    }

    //============================== Recompute k-eigenvalue
    double F_new = ComputeProduction();
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    //============================== Check convergence, reset book-keeping
    ScopedCopySTLvectors(groupset, phi_new_local, phi_prev_local);
    double k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    //============================== Print iteration summary
    if (options.verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info
          << chi_program_timer.GetTimeString() << " "
          << "  Iteration " << std::setw(5) << nit
          << "  k_eff " << std::setw(10) << k_eff
          << "  k_eff change " << std::setw(10) << k_eff_change
          << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi_log.Log(LOG_0) << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //============================== Initialize the precursor vector
  InitializePrecursors();

  double sweep_time = sweep_scheduler.GetAverageSweepTime();
  double source_time =
      chi_log.ProcessEvent(source_event_tag,
                           ChiLog::EventOperation::AVERAGE_DURATION);
  size_t num_angles = groupset.quadrature->abscissae.size();
  size_t num_unknowns = glob_node_count *
                        num_angles *
                        groupset.groups.size();
  chi_log.Log(LOG_0)
      << "\n";
  chi_log.Log(LOG_0)
      << "        Final k-eigenvalue    :        "
      << std::setprecision(6) << k_eff;
  chi_log.Log(LOG_0)
      << "        Set Src Time/sweep (s):        "
      << source_time;
  chi_log.Log(LOG_0)
      << "        Average sweep time (s):        "
      << sweep_time;
  chi_log.Log(LOG_0)
      << "        Sweep Time/Unknown (ns):       "
      << sweep_time * 1.0e9 * chi_mpi.process_count / static_cast<double>(num_unknowns);
  chi_log.Log(LOG_0)
      << "        Number of unknowns per sweep:  " << num_unknowns;
  chi_log.Log(LOG_0)
      << "\n\n";

  // GS_0 because we solve k-eigenvalue problems with 1 groupset right now.
  std::string sweep_log_file_name =
      std::string("GS_") + std::to_string(0) +
      std::string("_SweepLog_") + std::to_string(chi_mpi.location_id) +
      std::string(".log");
  groupset.PrintSweepInfoFile(sweep_scheduler.sweep_event_tag, sweep_log_file_name);

}
