#include <iostream>
#include "chi_runtime.h"

//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int main(int argc, char** argv)
{
  chi::run_time::Initialize(argc, argv);

  int error_code;
  if (chi::run_time::sim_option_interactive)
    error_code = chi::run_time::RunInteractive(argc, argv);
  else
    error_code = chi::run_time::RunBatch(argc, argv);

  chi::run_time::Finalize();

  return error_code;
}