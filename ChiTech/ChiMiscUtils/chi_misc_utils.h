#ifndef CHITECH_CHI_MISC_UTILS_H
#define CHITECH_CHI_MISC_UTILS_H

#include <cstddef>
#include <string>

/**Miscellaneous utilities. These utilities should have no dependencies.*/
namespace chi_misc_utils
{
  std::string PrintIterationProgress(size_t current_iteration,
                                     size_t total_num_iterations,
                                     unsigned int num_intvls = 10);
}//namespace chi_misc_utils

#endif //CHITECH_CHI_MISC_UTILS_H
