#include "chi_lua.h"

#include "chi_misc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_tests
{

int Test_chi_misc_utils(lua_State* L)
{
  bool passed = true;

  Chi::log.Log() << "Testing chi_misc_utils::PrintIterationProgress\n";

  const unsigned int I = 4;
  const size_t N=39;

  std::stringstream progress;
  for (size_t i=0; i<N; ++i)
  {
    progress << chi_misc_utils::PrintIterationProgress(i, N, I);
  }

  Chi::log.Log() << progress.str();

  return 0;
}

}//namespace chi_unit_tests