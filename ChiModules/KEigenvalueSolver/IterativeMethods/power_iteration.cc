#include "../k_eigenvalue_solver.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include <ChiTimer/chi_timer.h>

#include <iomanip>
#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepChunk SweepChunk;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

extern ChiTimer chi_program_timer;

using namespace LinearBoltzmann;

//###################################################################
/**Power iterative scheme for k-eigenvalue calculations.
 * Note that this routine currently only works when the problem 
 * is defined by a single groupset.
*/
void KEigenvalue::Solver::PowerIteration(LBSGroupset& groupset)
{
  chi_log.Log(LOG_0)
    << "\n\n";
  chi_log.Log(LOG_0)
    << "********** Solving k-eigenvalue problem with "
    << "the Power Method.\n\n";

  // ----- Setting up required sweep chunks
  auto sweep_chunk = SetSweepChunk(groupset);

  // ----- Set sweep scheduler
  MainSweepScheduler SweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                    &groupset.angle_agg);

  // ----- Tool the sweep chunk
  sweep_chunk->SetDestinationPhi(&phi_new_local);

  // ----- Set starting guess to a unit magnitude flux
  phi_prev_local.assign(phi_prev_local.size(),1.0);
  DisAssembleVectorLocalToLocal(groupset,phi_prev_local.data(),
                                         phi_old_local.data()); 

  // ----- Start outer k iterations
  double F_new = 0.0;        //Production source new (0.0 is not used)
  double F_prev = 1.0;       //Production source prev
  double k_eff_prev = k_eff;
  int nit = 0;               //number of iterations
  bool k_converged = false;

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
      SetKSource(groupset, APPLY_MATERIAL_SOURCE |
                           APPLY_SCATTER_SOURCE |
                           APPLY_FISSION_SOURCE);
      groupset.angle_agg.ZeroOutgoingDelayedPsi();
      phi_new_local.assign(phi_new_local.size(),0.0);
      SweepScheduler.Sweep(*sweep_chunk);

      // ----- Compute convergence parameters
      double pw_change = ComputePiecewiseChange(groupset);
      DisAssembleVectorLocalToLocal(groupset,phi_new_local.data(),
                                             phi_old_local.data());
      double rho = sqrt(pw_change/pw_change_prev);
      pw_change_prev = pw_change;
      nit += 1;

      if (si_nit==0) rho = 0.0;
      if (pw_change<std::max(groupset.residual_tolerance*rho,1.0e-10))
        si_converged = true;

      // ----- Print iteration information
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
    F_new = ComputeProduction();
    k_eff = F_new/F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    // ----- Compute convergence parameters and bump values
    double k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff; F_prev = F_new;
    DisAssembleVectorLocalToLocal(groupset,phi_new_local.data(),
                                           phi_prev_local.data());
                                          
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

  InitializePrecursors();
  
  double sweep_time = SweepScheduler.GetAverageSweepTime();
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
  groupset.PrintSweepInfoFile(SweepScheduler.sweep_event_tag,sweep_log_file_name);
  
}
