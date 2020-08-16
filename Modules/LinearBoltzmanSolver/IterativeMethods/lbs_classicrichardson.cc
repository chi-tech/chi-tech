#include "../lbs_linear_boltzman_solver.h"


#include "../../DiffusionSolver/Solver/diffusion_solver.h"

#include <ChiTimer/chi_timer.h>


#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <iomanip>

extern ChiTimer chi_program_timer;

//###################################################################
/**Solves a groupset using classic richardson.*/
void LinearBoltzman::Solver::ClassicRichardson(int group_set_num, SweepChunk* sweep_chunk,
    MainSweepScheduler & sweepScheduler, bool log_info /* = true*/)
{
  if (log_info)
  {
    chi_log.Log(LOG_0)
      << "\n\n";
    chi_log.Log(LOG_0)
      << "********** Solving groupset " << group_set_num
      << " with Classic-Richardson.\n\n";
  }

  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];
  int groupset_numgrps = groupset->groups.size();

  if (log_info)
    chi_log.Log(LOG_0)
      << "Quadrature number of angles: "
      << groupset->quadrature->abscissae.size() << "\n"
      << "Groups " << groupset->groups.front()->id << " "
      << groupset->groups.back()->id << "\n\n";

  //================================================== Tool the sweep chunk
  sweep_chunk->SetDestinationPhi(&phi_new_local);

  //================================================== Now start iterating
  double pw_change = 0.0;
  double pw_change_prev = 1.0;
  double rho = 0.0;
  bool converged = false;
  for (int k=0; k<groupset->max_iterations; k++)
  {
    SetSource(group_set_num,SourceFlags::USE_MATERIAL_SOURCE);

    groupset->angle_agg->ResetDelayedPsi();

    phi_new_local.assign(phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler.Sweep(sweep_chunk);

    if (groupset->apply_wgdsa)
    {
      AssembleWGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
      ((chi_diffusion::Solver*)groupset->wgdsa_solver)->ExecuteS(true,false);
      DisAssembleWGDSADeltaPhiVector(groupset, phi_new_local.data());
    }
    if (groupset->apply_tgdsa)
    {
      AssembleTGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
      ((chi_diffusion::Solver*)groupset->tgdsa_solver)->ExecuteS(true,false);
      DisAssembleTGDSADeltaPhiVector(groupset, phi_new_local.data());
    }

    pw_change = ComputePiecewiseChange(groupset);

    DisAssembleVectorLocalToLocal(groupset,phi_new_local.data(),
                                           phi_old_local.data());

    rho = sqrt(pw_change/pw_change_prev);
    pw_change_prev = pw_change;

    if (k==0) rho = 0.0;
    if (pw_change<std::max(groupset->residual_tolerance*rho,1.0e-10))
      converged = true;

    //======================================== Print iteration information
    std::string offset;
    if (groupset->apply_wgdsa || groupset->apply_tgdsa)
      offset = std::string("    ");

    std::stringstream iter_info;
    iter_info
      << chi_program_timer.GetTimeString() << " "
      << offset
      << "WGS groups ["
      << groupset->groups.front()->id
      << "-"
      << groupset->groups.back()->id
      << "]"
      << " Iteration " << std::setw(5) << k
      << " Point-wise change " << std::setw(14) << pw_change;

    if (converged)
      iter_info << " CONVERGED\n";

    if (log_info)
      chi_log.Log(LOG_0) << iter_info.str();

    if (converged) break;

    if (options.write_restart_data)
    {
      if ((chi_program_timer.GetTime()/60000.0) >
          last_restart_write+options.write_restart_interval)
      {
        last_restart_write = chi_program_timer.GetTime()/60000.0;
        WriteRestartData(options.write_restart_folder_name,
                         options.write_restart_file_base);
      }
    }
  }


  double sweep_time = sweepScheduler.GetAverageSweepTime();
  double source_time=
    chi_log.ProcessEvent(source_event_tag,
                         ChiLog::EventOperation::AVERAGE_DURATION);
  size_t num_angles = groupset->quadrature->abscissae.size();
  long int num_unknowns = (long int)glob_dof_count*
                          (long int)num_angles*
                          (long int)groupset->groups.size();

  if (log_info)
  {
    chi_log.Log(LOG_0)
      << "\n\n";
    chi_log.Log(LOG_0)
      << "        Set Src Time/sweep (s):        "
      << source_time;
    chi_log.Log(LOG_0)
      << "        Average sweep time (s):        "
      << sweep_time;
    chi_log.Log(LOG_0)
      << "        Sweep Time/Unknown (ns):       "
      << sweep_time*1.0e9*chi_mpi.process_count/num_unknowns;
    chi_log.Log(LOG_0)
      << "        Number of unknowns per sweep:  " << num_unknowns;
    chi_log.Log(LOG_0)
      << "\n\n";
  }

  std::string sweep_log_file_name =
    std::string("GS_") + std::to_string(group_set_num) +
    std::string("_SweepLog_") + std::to_string(chi_mpi.location_id) +
    std::string(".log");
  groupset->PrintSweepInfoFile(sweepScheduler.sweep_event_tag,sweep_log_file_name);
}
