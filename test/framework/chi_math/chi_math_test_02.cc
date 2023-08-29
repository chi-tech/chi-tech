#include "math/ParallelVector/parallel_vector.h"
#include "math/ParallelVector/ghosted_parallel_vector.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

namespace chi_unit_tests
{

chi::ParameterBlock
chi_math_Test02_ParallelVector(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_Test02_ParallelVector,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_math_Test02_ParallelVector);

chi::ParameterBlock
chi_math_Test02_ParallelVector(const chi::InputParameters& params)
{
  using namespace chi_math;

  ChiLogicalErrorIf(Chi::mpi.process_count != 2, "Requires 2 processors");

  Chi::log.Log() << "Testing chi_math::ParallelVector" << std::endl;

  ParallelVector vec(5, 10, Chi::mpi.comm);

  if (Chi::mpi.location_id == 0) vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    vec.SetValue(0, 1.0, VecOpType::SET_VALUE);
  vec.Assemble();

  Chi::log.LogAll() << "vec after assembly: " << vec.PrintStr() << std::endl;

  Chi::log.Log() << "Testing chi_math::GhostedParallelVector" << std::endl;

  const int64_t ghost_id = Chi::mpi.location_id == 0 ? 5 : 4;
  GhostedParallelVector ghost_vec(5, 10, {ghost_id}, Chi::mpi.comm);

  Chi::log.LogAll() << "Number of Ghosts: " << ghost_vec.NumGhosts()
                    << std::endl;

  if (Chi::mpi.location_id == 0)
    ghost_vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    ghost_vec.SetValue(4, 1.0, VecOpType::SET_VALUE);
  ghost_vec.Assemble();
  ghost_vec.CommunicateGhostEntries();

  Chi::log.LogAll() << "Ghost vec after communicate: " << ghost_vec.PrintStr()
                    << std::endl;

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
