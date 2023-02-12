#include "lbts_transient_solver.h"

#include "B_LBSSteadyState/IterativeOperations/sweep_wgs_context.h"
#include "B_LBSSteadyState/IterativeMethods/wgs_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Transient solver execute routine.*/
void lbs::TransientSolver::Execute()
{
  chi::log.Log() << "Executing " << TextName() << ".";

  const int max_num_steps = transient_options.max_time_steps;
  const double max_time = transient_options.t_final;
  int step_number = 0;
  while (((max_num_steps > 0 and step_number < max_num_steps) or
         (max_num_steps < 0)) and (time < max_time))
  {
    Step();

    PostStepCallBackFunction();

    if (not transient_options.inhibit_advance)
    {
      AdvanceTimeValues(); //new copied to prev + time+=dt
      ++step_number;
      transient_options.inhibit_advance = false;
    }
  }

  UpdateFieldFunctions();

  chi::log.Log() << "Done Executing " << TextName() << ".";
}

//###################################################################
/**Transient solver timestep routine.*/
void lbs::TransientSolver::Step()
{
  if (transient_options.verbosity_level >= 2)
    chi::log.Log() << TextName() << " Stepping with dt " << dt;

  phi_old_local_ = phi_prev_local;

  for (auto& groupset : groupsets_)
  {
    //======================================== Converge the scattering source
    //                                         with a fixed fission source
    //                                         and temporal source
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    auto sweep_chunk = SetTransientSweepChunk(groupset);

    auto sweep_wgs_context_ptr =
      std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
        *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,  //lhs_scope
        APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
        APPLY_AGS_FISSION_SOURCES,                              //rhs_scope
        true/*with_delayed_psi*/,
        options_.verbose_inner_iterations,
        sweep_chunk);

    WGSLinearSolver<Mat,Vec,KSP> solver(sweep_wgs_context_ptr);
    solver.Setup();
    solver.Solve();

    MPI_Barrier(MPI_COMM_WORLD);
  }

  //======================================== Compute t^{n+1} value
  {
    const auto& BackwardEuler = chi_math::SteppingMethod::BACKWARD_EULER;
    const auto& CrankNicolson = chi_math::SteppingMethod::CRANK_NICHOLSON;

    double theta;
    if      (method == BackwardEuler) theta = 1.0;
    else if (method == CrankNicolson) theta = 0.5;
    else                              theta = 0.7;
    const double inv_theta = 1.0/theta;

    auto& phi = phi_new_local_;
    const auto& phi_prev = phi_prev_local;
    for (size_t i = 0; i < phi.size(); ++i)
      phi[i] = inv_theta*(phi[i] + (theta-1.0) * phi_prev[i]);

    if (options_.use_precursors)
      StepPrecursors();
  }

  const double FR_new = ComputeFissionProduction(phi_new_local_);

  //============================================= Print end of timestep
  if (transient_options.verbosity_level >= 1)
  {
    char buff[200];
    snprintf(buff,200, " dt=%.1e time=%10.4g FR=%12.6g", dt, time + dt, FR_new);
    chi::log.Log() << TextName() << buff;
  }

  UpdateFieldFunctions();
}


//###################################################################
/**Advances time values.*/
void lbs::TransientSolver::AdvanceTimeValues()
{
  time += dt;
  phi_prev_local = phi_new_local_;
  psi_prev_local = psi_new_local_;
  if (options_.use_precursors)
    precursor_prev_local = precursor_new_local_;
}