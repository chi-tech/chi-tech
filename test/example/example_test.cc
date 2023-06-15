#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_tests
{

chi_objects::ParameterBlock ExampleTest(const chi_objects::InputParameters&);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/ExampleTest,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/ExampleTest);

chi_objects::ParameterBlock ExampleTest(const chi_objects::InputParameters&)
{
  chi::log.Log() << "This is an example test";

  return chi_objects::ParameterBlock();
}

} // namespace chi_unit_tests