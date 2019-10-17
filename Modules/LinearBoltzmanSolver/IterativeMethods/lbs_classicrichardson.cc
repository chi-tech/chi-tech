#include "../lbs_linear_boltzman_solver.h"

#include "ChiMesh/SweepUtilities/chi_sweepscheduler.h"
#include "../SweepChunks/lbs_sweepchunk_pwl_polyhedron.h"
#include "../../DiffusionSolver/Solver/diffusion_solver.h"

#include <ChiTimer/chi_timer.h>


#include <chi_log.h>

extern ChiLog chi_log;

extern double chi_global_timings[20];

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;
typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;

extern ChiTimer chi_program_timer;

//###################################################################
/**Solves a groupset using classic richardson.*/
void LinearBoltzman::Solver::ClassicRichardson(int group_set_num)
{
  chi_log.Log(LOG_0)
    << "\n\n";
  chi_log.Log(LOG_0)
    << "********** Solving groupset " << group_set_num
    << " with Classic-Richardson.\n\n";

  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];
  int groupset_numgrps = groupset->groups.size();
  chi_log.Log(LOG_0)
    << "Quadrature number of angles: "
    << groupset->quadrature->abscissae.size() << "\n"
    << "Number of azimuthal angles : "
    << groupset->quadrature->azimu_ang.size() << "\n"
    << "Number of polar angles     : "
    << groupset->quadrature->polar_ang.size() << "\n"
    << "Groups " << groupset->groups.front()->id << " "
    << groupset->groups.back()->id << "\n\n";

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = SetSweepChunk(group_set_num);

  //================================================== Set sweep scheduler
  MainSweepScheduler sweepScheduler(DEPTH_OF_GRAPH,
                                    groupset->angle_agg);

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

    groupset->latest_convergence_metric = std::min(pw_change, 1.0);
    ConvergeCycles(sweepScheduler,sweep_chunk,groupset);

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

    chi_log.Log(LOG_0) << iter_info.str();

    if (converged) break;
  }


  double sweep_time = chi_global_timings[16]/chi_global_timings[17];
  double source_time= chi_global_timings[18]/chi_global_timings[19];
  size_t num_angles = groupset->quadrature->abscissae.size();
  long int num_unknowns = (long int)glob_dof_count*
                          (long int)num_angles*
                          (long int)groupset->groups.size();
  chi_log.Log(LOG_0)
    << "\n\n";
  chi_log.Log(LOG_0)
    << "        Set Src Time/Unknown (ns):    "
    << source_time*1.0e9*chi_mpi.process_count/num_unknowns;
  chi_log.Log(LOG_0)
    << "        Sweep Time/Unknown (ns):       "
    << sweep_time*1.0e9*chi_mpi.process_count/num_unknowns;
  chi_log.Log(LOG_0)
    << "        Number of unknowns per sweep:  " << num_unknowns;
  chi_log.Log(LOG_0)
    << "\n\n";
}