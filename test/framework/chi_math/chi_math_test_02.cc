#include "math/ParallelVector/parallel_vector.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#include <unistd.h>


namespace chi_unit_tests
{

chi::ParameterBlock
chi_math_Test02(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_Test02,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_math_Test02);

chi::ParameterBlock
chi_math_Test02(const chi::InputParameters& params)
{
  using namespace chi_math;

  Chi::log.Log() << "GOLD_BEGIN";

  Chi::log.Log() << "Testing chi_math::ParallelVector" << std::endl;

  ParallelVector vec(5, 10, Chi::mpi.comm);

  if (Chi::mpi.location_id == 0)
    vec.SetValue(5, 2.0, OperationType::SET_VALUE);
  else
    vec.SetValue(0, 1.0, OperationType::SET_VALUE);
  vec.Assemble();

  Chi::mpi.Barrier();
  Chi::log.LogAll() << vec.PrintStr() << std::endl;
  Chi::mpi.Barrier();
  usleep(1'000);

  Chi::log.Log()
      << "Testing chi_math::GhostedParallelVector" << std::endl;

  const int64_t ghost_id = Chi::mpi.location_id == 0 ? 5 : 4;
  GhostedParallelVector ghost_vec(5, 10, {ghost_id}, Chi::mpi.comm);

  Chi::mpi.Barrier();
  Chi::log.LogAll()
      << "Number of Ghosts: " << ghost_vec.NumGhosts() << std::endl;
  Chi::mpi.Barrier();
  usleep(1'000);

  if (Chi::mpi.location_id == 0)
    ghost_vec.SetValue(5, 2.0, OperationType::SET_VALUE);
  else
    ghost_vec.SetValue(4, 1.0, OperationType::SET_VALUE);
  ghost_vec.Assemble();
  ghost_vec.CommunicateGhostEntries();

  Chi::mpi.Barrier();
  Chi::log.LogAll() << ghost_vec.PrintStr() << std::endl;
  Chi::mpi.Barrier();
  usleep(1'000);

  Chi::log.Log() << "GOLD_END";
  return chi::ParameterBlock();
}

}
