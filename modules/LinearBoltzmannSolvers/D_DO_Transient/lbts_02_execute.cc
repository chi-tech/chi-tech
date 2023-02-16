#include "lbts_transient_solver.h"

#include "chi_log.h"

//###################################################################
/**Transient solver execute routine.*/
void lbs::DiscOrdTransientSolver::Execute()
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
