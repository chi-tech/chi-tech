#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

#include "../Source/ResidualSource/mc_rmc_source.h"

extern ChiLog chi_log;
extern ChiTimer chi_program_timer;
typedef unsigned long long TULL;


extern ChiMath chi_math_handler;

//#########################################################
/**Executes the solver*/
void chi_montecarlon::Solver::Execute()
{
  chi_log.Log(LOG_0) << "Executing Montecarlo solver";

  int block_lengths[] = {3,3,1,2,2};
  MPI_Aint block_disp[] = {
    offsetof(Particle,pos),
    offsetof(Particle,dir),
    offsetof(Particle,w),
    offsetof(Particle,egrp),
    offsetof(Particle,alive)};

  MPI_Datatype types[] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_INT,
    MPIU_BOOL
  };

  MPI_Datatype prtcl_data_type;

  MPI_Type_create_struct(5,block_lengths,block_disp,types,&prtcl_data_type);

  chi_montecarlon::Source* src = sources.back();

  nps_global = 0;
  TULL nps_last = 0;
  double start_time = chi_program_timer.GetTime()/1000.0;
  double time = 0.0;
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

      while (prtcl.alive and !prtcl.banked)
        Raytrace(&prtcl);

    }//for pi in batch

    MPI_Barrier(MPI_COMM_WORLD);



    RendesvouzTallies();
    if (make_pwld)
      RendesvouzPWLTallies();
    ComputeRelativeStdDev();



    time = chi_program_timer.GetTime()/1000.0;
    particle_rate = (nps_global)*3600.0e-6/(time-start_time);
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString()
      << " TFC-rendesvouz: # of particles ="
      << std::setw(14)
      << nps_global
      << " avg-rate = "
      << std::setw(6) << std::setprecision(4)
      << particle_rate << " M/hr"
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
  NormalizeTallies();
  if (make_pwld)
    NormalizePWLTallies();

  //==================== Computing pwl transformations
  size_t num_cells = grid->local_cell_glob_indices.size();
  if (make_pwld)
  {
    for (size_t lc=0; lc<num_cells; lc++)
    {
      int cell_g_index = grid->local_cell_glob_indices[lc];
      int map = local_cell_pwl_dof_array_address[lc];

      auto cell_pwl_view = pwl_discretization->MapFeView(cell_g_index);

      MatDbl A(cell_pwl_view->IntV_shapeI_shapeJ);
      MatDbl Ainv = chi_math_handler.Inverse(A);
      VecDbl b(cell_pwl_view->dofs,0.0);

      for (int g=0; g<num_grps; g++)
      {
        for (int dof=0; dof<cell_pwl_view->dofs; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          b[dof] = phi_pwl_global[ir]*cell_pwl_view->IntV_shapeI[dof];
        }
        VecDbl x = chi_math_handler.MatMul(Ainv,b);
        for (int dof=0; dof<cell_pwl_view->dofs; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          phi_pwl_global[ir] = x[dof];
        }
      }
    }//for local cell lc
  }//if make_pwld


  //Print group 0

//  for (int lc=0; lc<num_cells; lc++)
//  {
//    std::stringstream outstr;
//    outstr
//      << "Cell " << lc
//      << " phi=" << phi_global[lc*num_grps]
//      << " std=" << phi_local_relsigma[lc*num_grps+0]
//      << " abs=" << phi_local_relsigma[lc*num_grps+0]*phi_global[lc*num_grps];
//
//    if (make_pwld)
//    {
//      int map0 = local_cell_pwl_dof_array_address[lc];
//      int map1 = map0 + 1*num_grps*num_moms + num_grps*0 + 0;
//
//      outstr
//        << " dof0g0m0=" << phi_pwl_global[map0]
//        << " dof1g0m0=" << phi_pwl_global[map1];
//    }
//
//
//    chi_log.Log(LOG_0) << outstr.str();
//  }


  auto rsrc = (chi_montecarlon::ResidualSource*)sources[0];
  std::cout << "Particles sampled intr = " << rsrc->particles_C << std::endl;
  std::cout << "Particles sampled left = " << rsrc->particles_L << std::endl;
  std::cout << "Particles sampled rite = " << rsrc->particles_R << std::endl;
  std::cout << "Weights sampled left = " << rsrc->weights_L << std::endl;
  std::cout << "Weights sampled rite = " << rsrc->weights_R << std::endl;

  chi_log.Log(LOG_0) << "Done executing Montecarlo solver";
}