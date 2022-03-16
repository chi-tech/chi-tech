#include "unit_tests.h"

#include "chi_misc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;


bool chi_unit_tests::Test_chi_misc_utils(bool verbose)
{
  bool passed = true;
  std::stringstream output;

  output << "Testing chi_misc_utils::PrintIterationProgress\n";

  const unsigned int I = 4;
  const size_t N=39;

  std::stringstream progress;
  for (size_t i=0; i<N; ++i)
  {
    progress << chi_misc_utils::PrintIterationProgress(i, N, I);
  }

  if (progress.str() != "  25.00  50.00  75.00 100.00")
  {
    passed = false;
    output << std::string("chi_misc_utils::PrintIterationProgress(.,39,4) ... Failed\n"
                          " Expected:\n"
                          "  25.00  50.00  75.00 100.00\n"
                          "Got:\n") + progress.str();
  }
  else
    output << std::string("chi_misc_utils::PrintIterationProgress(.,39,4) ... Passed\n");

  if (verbose)
    chi_log.Log() << output.str();

  return passed;
}