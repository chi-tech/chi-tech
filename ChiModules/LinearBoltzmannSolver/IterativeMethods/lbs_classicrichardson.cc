#include "../lbs_linear_boltzmann_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

#include <iomanip>

#define sc_double static_cast<double>

//###################################################################
/**Solves a groupset using classic richardson.*/
bool lbs::SteadySolver::
ClassicRichardson(LBSGroupset& groupset,
                  MainSweepScheduler& sweep_scheduler,
                  SourceFlags source_flags,
                  bool log_info /* = true*/)
{
  constexpr bool WITH_DELAYED_PSI = true;
  if (log_info)
  {
    chi::log.Log() << "\n\n";
    chi::log.Log() << "********** Solving groupset" << groupset.id
                       << " with Classic-Richardson.\n\n";
    chi::log.Log()
      << "Quadrature number of angles: "
      << groupset.quadrature->abscissae.size() << "\n"
      << "Groups " << groupset.groups.front().id << " "
      << groupset.groups.back().id << "\n\n";
  }

  const auto num_delayed_psi_info = groupset.angle_agg.GetNumDelayedAngularDOFs();
  const size_t num_angles = groupset.quadrature->abscissae.size();
  const size_t num_psi_global = glob_node_count *
                                num_angles *
                                groupset.groups.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info)
  {
    chi::log.Log()
      << "Total number of angular unknowns: "
      << num_psi_global
      << "\n"
      << "Number of lagged angular unknowns: "
      << num_delayed_psi_globl << "("
      << sc_double(num_delayed_psi_globl) / sc_double(num_psi_global)
      << "%)";
  }

  std::vector<double> init_q_moments_local = q_moments_local;

  //================================================== Sweep chunk settings
  auto& sweep_chunk = sweep_scheduler.GetSweepChunk();
  bool use_surface_source_flag = (source_flags & APPLY_MATERIAL_SOURCE) and
                                 (not options.use_src_moments);
  sweep_chunk.SetSurfaceSourceActiveFlag(use_surface_source_flag);
  sweep_chunk.ZeroIncomingDelayedPsi();

  //================================================== Now start iterating
  double pw_change_prev = 1.0;
  bool converged = false;
  for (int k = 0; k < groupset.max_iterations; ++k)
  {
    q_moments_local = init_q_moments_local;
    SetSource(groupset, q_moments_local, source_flags);

    sweep_chunk.ZeroFluxDataStructures();
    sweep_scheduler.Sweep();

    if (groupset.apply_wgdsa)
      ExecuteWGDSA(groupset,phi_old_local,phi_new_local);

    if (groupset.apply_tgdsa)
      ExecuteTGDSA(groupset,phi_old_local,phi_new_local);

    double pw_change = ComputePiecewiseChange(groupset);

    ScopedCopySTLvectors(groupset,phi_new_local,phi_old_local, WITH_DELAYED_PSI);

    double rho = sqrt(pw_change / pw_change_prev);
    pw_change_prev = pw_change;

    if (k==0) rho = 0.0;

    if (pw_change<std::max(groupset.residual_tolerance*(1.0-rho),1.0e-10))
      converged = true;

    //======================================== Print iteration information
    {
      std::string offset;
      if (groupset.apply_wgdsa || groupset.apply_tgdsa)
        offset = std::string("    ");

    std::stringstream iter_info;
    iter_info
      << chi::program_timer.GetTimeString() << " "
      << offset
      << "WGS groups ["
      << groupset.groups.front().id
      << "-"
      << groupset.groups.back().id
      << "]"
      << " Iteration " << std::setw(5) << k
      << " Point-wise change " << std::setw(14) << pw_change
      <<" Spectral Radius Estimate " << std::setw(10) << rho;

      if (converged)
        iter_info << " CONVERGED\n";

      if (log_info)
        chi::log.Log() << iter_info.str();

      if (converged) break;

      if (options.write_restart_data)
      {
        if ((chi::program_timer.GetTime()/60000.0) >
            last_restart_write+options.write_restart_interval)
        {
          last_restart_write = chi::program_timer.GetTime()/60000.0;
          WriteRestartData(options.write_restart_folder_name,
                           options.write_restart_file_base);
        }
      }//if write restart data
    }//print iterative info

  }

  //============================================= Print solution info
  {
    double sweep_time = sweep_scheduler.GetAverageSweepTime();
    double source_time=
      chi::log.ProcessEvent(source_event_tag,
                           chi_objects::ChiLog::EventOperation::AVERAGE_DURATION);

    if (log_info)
    {
      chi::log.Log()
        << "\n\n";
      chi::log.Log()
        << "        Set Src Time/sweep (s):        "
        << source_time;
      chi::log.Log()
        << "        Average sweep time (s):        "
        << sweep_time;
      chi::log.Log()
        << "        Sweep Time/Unknown (ns):       "
        << sweep_time*1.0e9*chi::mpi.process_count/
            sc_double(num_psi_global);
      chi::log.Log()
        << "        Number of unknowns per sweep:  " << num_psi_global;
      chi::log.Log()
        << "\n\n";

      std::string sweep_log_file_name =
          std::string("GS_") + std::to_string(groupset.id) +
          std::string("_SweepLog_") + std::to_string(chi::mpi.location_id) +
          std::string(".log");
      groupset.PrintSweepInfoFile(sweep_scheduler.sweep_event_tag, sweep_log_file_name);
    }
  }//print solution info

  return converged;
}
