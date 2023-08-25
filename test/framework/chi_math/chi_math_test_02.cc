#include "math/parallel_vector.h"
#include "math/ghosted_parallel_vector.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"


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
  Chi::log.Log() << "GOLD_BEGIN";

  Chi::mpi.Barrier();
  if (Chi::mpi.location_id == 0) Chi::mpi.Barrier();
  Chi::log.LogAll() << "Testing chi_math::ParallelVector" << std::endl;
  Chi::mpi.Barrier();

  chi_math::ParallelVector vec(5, 10, Chi::mpi.comm);

  if (Chi::mpi.location_id == 0)
    vec.AddValues({5, 1}, {2.0, 2.0});
  else
    vec.AddValues({0, 6}, {1.0, 1.0});
  vec.Assemble();

  Chi::mpi.Barrier();
  if (Chi::mpi.location_id == 0) Chi::mpi.Barrier();
  Chi::log.LogAll() << vec.PrintStr() << std::endl;
  Chi::mpi.Barrier();

  Chi::mpi.Barrier();
  if (Chi::mpi.location_id == 0) Chi::mpi.Barrier();
  Chi::log.LogAll()
      << "Testing chi_math::GhostedParallelVector" << std::endl;
  Chi::mpi.Barrier();

  const int64_t ghost_id = Chi::mpi.location_id == 0 ? 5 : 4;
  chi_math::GhostedParallelVector ghost_vec(5, 10, {ghost_id}, Chi::mpi.comm);

  Chi::mpi.Barrier();
  if (Chi::mpi.location_id == 0) Chi::mpi.Barrier();
  Chi::log.LogAll()
      << "Number of Ghosts: " << ghost_vec.NumGhosts() << std::endl;
  Chi::mpi.Barrier();

  if (Chi::mpi.location_id == 0)
    ghost_vec.AddValue(5, 2.0);
  else
    ghost_vec.AddValue(4, 1.0);
  ghost_vec.Assemble();

  Chi::mpi.Barrier();
  if (Chi::mpi.location_id == 0) Chi::mpi.Barrier();
  Chi::log.LogAll() << ghost_vec.PrintStr() << std::endl;
  Chi::mpi.Barrier();

  Chi::log.Log() << "GOLD_END";
  return chi::ParameterBlock();
}

}
