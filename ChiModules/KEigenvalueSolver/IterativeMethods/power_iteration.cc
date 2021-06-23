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
    << "\n\n";
  chi_log.Log(LOG_0)
    << "********** Solving k-eigenvalue problem with "
    << "the Power Method.\n\n";

  LBSGroupset&  groupset = group_sets[0];

  groupset.angle_agg.ZeroIncomingDelayedPsi();

  //======================================== Setup sweep chunk
  auto sweep_chunk = SetSweepChunk(groupset);
  MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                     groupset.angle_agg,
                                     *sweep_chunk);

  //======================================== Tool the sweep chunk
  sweep_scheduler.sweep_chunk.SetDestinationPhi(phi_new_local);

  //======================================== Initial guess
  phi_prev_local.assign(phi_prev_local.size(),1.0);
  ScopedCopySTLvectors(groupset, phi_prev_local, phi_old_local);

  //======================================== Start power iterations
  double F_prev     = 1.0;    //production source prev
  double k_eff_prev = k_eff;
  int nit           = 0;      //number of iterations
  bool k_converged  = false;

  while (nit < options.max_iterations)
  {
    chi_log.Log(LOG_0VERBOSE_2)
      << "\n********** Starting source iterations";

    // ----- Start inner source iterations
    double pw_change_prev = 1.0;
    bool si_converged = false;
    for (int si_nit=0; si_nit<groupset.max_iterations; si_nit++)
    {
      // ----- Set source and sweep
      q_moments_local.assign(q_moments_local.size(), 0.0);
      SetKSource(groupset, q_moments_local,
                 APPLY_AGS_SCATTER_SOURCE | APPLY_WGS_SCATTER_SOURCE |
                 APPLY_AGS_FISSION_SOURCE | APPLY_WGS_FISSION_SOURCE);

      groupset.ZeroAngularFluxDataStructures();
      phi_new_local.assign(phi_new_local.size(),0.0);
      sweep_scheduler.Sweep();

      //======================================== Compute convergence params
      double pw_change = ComputePiecewiseChange(groupset);
      ScopedCopySTLvectors(groupset, phi_new_local, phi_old_local);
      double rho = sqrt(pw_change/pw_change_prev);
      pw_change_prev = pw_change;
      nit += 1;

      if (si_nit==0) rho = 0.0;
      if (pw_change<std::max(groupset.residual_tolerance*rho,1.0e-10))
        si_converged = true;

      //============================== Print iteration information
      std::string offset = "    ";
      std::stringstream si_iter_info;
      si_iter_info
        << chi_program_timer.GetTimeString() << " "
        << offset
        << "WGS groups ["
        << groupset.groups.front().id
        << "-"
        << groupset.groups.back().id
        << "]"
        << " Source Iteration " << std::setw(5) << si_nit
        << " Point-wise change " << std::setw(14) << pw_change;

      if (si_converged)
        si_iter_info << " CONVERGED\n";
      chi_log.Log(LOG_0VERBOSE_1) << si_iter_info.str();

      if (si_converged) {
        chi_log.Log(LOG_0VERBOSE_1)
            << "\nSource iterations converged in "
            << si_nit << " iteratrions.\n";
        break;
      }
    }//for source iterations

    if (!si_converged) {
      chi_log.Log(LOG_ALLVERBOSE_1)
          << "\n!!!WARNING!!! "
             "Source iterations did not converge.\n";
    }

    // ----- Recompute eigenvalue
    double F_new = ComputeProduction();
    k_eff = F_new/F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    // ----- Compute convergence parameters and bump values
    double k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff; F_prev = F_new;
    ScopedCopySTLvectors(groupset, phi_new_local, phi_prev_local);
                                          
    if (k_eff_change<std::max(options.tolerance, 1.0e-12))
      k_converged = true;    

    // ----- Print iteration summary
    std::stringstream k_iter_info;
    k_iter_info
      << chi_program_timer.GetTimeString() << " "
      << "  Iteration " << std::setw(5) << nit
      << "  k_eff " << std::setw(10) << k_eff
      << "  k_eff change " << std::setw(10) << k_eff_change
      << "  reactivity " << std::setw(10) << reactivity * 1e5;
    if (k_converged) {
      k_iter_info << " CONVERGED\n";
    }
    chi_log.Log(LOG_0VERBOSE_1) << k_iter_info.str();

    if (k_converged) break;
  }//for k iterations

  if (options.use_precursors)
    InitializePrecursors();
  
  double sweep_time = sweep_scheduler.GetAverageSweepTime();
  double source_time=
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
    << sweep_time*1.0e9*chi_mpi.process_count/static_cast<double>(num_unknowns);
  chi_log.Log(LOG_0)
    << "        Number of unknowns per sweep:  " << num_unknowns;
  chi_log.Log(LOG_0)
    << "\n\n";

  // GS_0 because we solve k-eigenvalue problems with 1 groupset right now.
  std::string sweep_log_file_name =
    std::string("GS_") + std::to_string(0) +
    std::string("_SweepLog_") + std::to_string(chi_mpi.location_id) +
    std::string(".log");
  groupset.PrintSweepInfoFile(sweep_scheduler.sweep_event_tag,sweep_log_file_name);
  
}
