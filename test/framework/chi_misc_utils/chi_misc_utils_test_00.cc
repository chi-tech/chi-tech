#include "chi_misc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiConsole/chi_console.h"

namespace chi_unit_tests
{

chi_objects::ParameterBlock
chi_misc_utils_Test00(const chi_objects::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_misc_utils_Test00,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_misc_utils_Test00);

chi_objects::ParameterBlock
chi_misc_utils_Test00(const chi_objects::InputParameters&)
{
  chi::log.Log() << "Testing chi_misc_utils::PrintIterationProgress\n";

  const unsigned int I = 4;
  const size_t N = 39;

  std::stringstream progress;
  for (size_t i = 0; i < N; ++i)
  {
    progress << chi_misc_utils::PrintIterationProgress(i, N, I);
  }

  chi::log.Log() << progress.str();

  return chi_objects::ParameterBlock();
}

} // namespace chi_unit_tests