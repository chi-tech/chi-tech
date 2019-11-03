#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

#include "../Source/ResidualSource/mc_rmc_source.h"

#include <PiecewiseLinear/CellViews/pwl_slab.h>

extern ChiLog chi_log;
extern ChiTimer chi_program_timer;
typedef unsigned long long TULL;

#include<typeinfo>

extern ChiMath chi_math_handler;

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
      if (std::isnan(prtcl.dir.x))
      {
        chi_log.Log(LOG_ALLERROR)
          << "Particle dir corrupt.";
        exit(EXIT_FAILURE);
      }

      while (prtcl.alive)
      {
        Raytrace(&prtcl);
//        chi_log.Log(LOG_0)
//          << "Ray particle " << prtcl.pos.PrintS()
//          << " " << prtcl.dir << " " << prtcl.cur_cell_ind;
//        usleep(100000);
      }


      ComputeTallySqr();
      ComputePWLTallySqr();
    }//for pi in batch




    RendesvouzTallies();
    RendesvouzPWLTallies();
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

  for (auto& tally_value : phi_pwl_global)
    tally_value *= tally_multipl_factor/nps_global;

  //==================== Computing pwl transformations
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_g_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_g_index];
    int map = local_cell_pwl_dof_array_address[lc];

    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto cell_pwl_view =
        static_cast<SlabFEView*>(pwl_discretization->MapFeView(cell_g_index));

      MatDbl A(cell_pwl_view->IntV_shapeI_shapeJ);
      MatDbl Ainv = chi_math_handler.Inverse(A);
      VecDbl b(2,0.0);

      for (int g=0; g<num_grps; g++)
      {
        for (int dof=0; dof<2; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          b[dof] = phi_pwl_global[ir]*cell_pwl_view->IntV_shapeI[dof];
        }
        VecDbl x = chi_math_handler.MatMul(Ainv,b);
        for (int dof=0; dof<2; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          phi_pwl_global[ir] = x[dof];
        }
      }


    }
  }

  //Print group 0

  for (int lc=0; lc<num_cells; lc++)
  {
    int map0 = local_cell_pwl_dof_array_address[lc];
    int map1 = map0 + 1*num_grps*num_moms + num_grps*0 + 0;;
    chi_log.Log(LOG_0)
      << "Cell " << lc
      << " phi=" << phi_global[lc*num_grps]
      << " std=" << phi_local_relsigma[lc*num_grps+0]
      << " abs=" << phi_local_relsigma[lc*num_grps+0]*phi_global[lc*num_grps]
      << " dof0g0m0=" << phi_pwl_global[map0]
      << " dof1g0m0=" << phi_pwl_global[map1];
  }


  chi_montecarlon::ResidualSource* rsrc = (chi_montecarlon::ResidualSource*)sources[0];
  std::cout << "Particles sampled intr = " << rsrc->particles_C << std::endl;
  std::cout << "Particles sampled left = " << rsrc->particles_L << std::endl;
  std::cout << "Particles sampled rite = " << rsrc->particles_R << std::endl;
  std::cout << "Weights sampled left = " << rsrc->weights_L << std::endl;
  std::cout << "Weights sampled rite = " << rsrc->weights_R << std::endl;

  chi_log.Log(LOG_0) << "Done executing Montecarlo solver";
}