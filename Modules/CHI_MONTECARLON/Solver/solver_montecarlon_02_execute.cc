#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiLog chi_log;
extern ChiTimer chi_program_timer;
typedef unsigned long long TULL;

//#########################################################
/**Executes the solver*/
void chi_montecarlon::Solver::Execute()
{
  chi_log.Log(LOG_0) << "Executing Montecarlo solver";

  chi_montecarlon::Source* src = sources.back();

  nps_global = 0;
  TULL nps_last = 0;
  double start_time = chi_program_timer.GetTime()/1000.0;
  double time = 0.0;
  double last_time = 0.0;
  double particle_rate = 0.0;
  for (int b=0; b<batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    nps_last = 0;
    current_batch = b;
    for (TULL pi=0; pi<batch_sizes_per_loc[b]; pi++)
    {
      nps++;
      nps_last++;
      chi_montecarlon::Particle prtcl = src->CreateParticle(&rng0);

//      chi_log.Log(LOG_0)
//        << "Src particle " << prtcl.pos.PrintS()
//        << " " << prtcl.dir << " " << prtcl.cur_cell_ind;
//      usleep(100000);

      while (prtcl.alive)
      {
        Raytrace(&prtcl);
//        chi_log.Log(LOG_0)
//          << "Ray particle " << prtcl.pos.PrintS()
//          << " " << prtcl.dir << " " << prtcl.cur_cell_ind;
//        usleep(100000);
      }


      ComputeTallySqr();
    }//for pi in batch




    RendesvouzTallies();
    ComputeRelativeStdDev();



    time = chi_program_timer.GetTime()/1000.0;
    particle_rate = (nps_global)*60.0e-6/(time-start_time);
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString()
      << " TFC-rendesvouz: # of particles ="
      << std::setw(14)
      << nps_global
      << " avg-rate = "
      << std::setw(6) << std::setprecision(4)
      << particle_rate << " M/min"
      << " Max Rel.Sigma = "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error
      << " "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error2
      << " "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error3;

  }

  //Normalize tallies
  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    for (int g=0; g<num_grps; g++)
    {
      int ir = lc*num_grps + g;
      phi_global[ir] = phi_global[ir]*tally_multipl_factor/nps_global;
    }
  }

  //Print group 0

  for (int lc=0; lc<num_cells; lc++)
  {
    chi_log.Log(LOG_0)
      << "Cell " << lc
      << " phi=" << phi_global[lc*num_grps]
      << " std=" << phi_local_relsigma[lc*num_grps+0];
  }


  chi_log.Log(LOG_0) << "Done executing Montecarlo solver";
}