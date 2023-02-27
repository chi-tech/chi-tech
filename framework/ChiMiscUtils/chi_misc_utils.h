#ifndef CHITECH_CHI_MISC_UTILS_H
#define CHITECH_CHI_MISC_UTILS_H

#include <cstddef>
#include <string>
#include <algorithm>

/**Miscellaneous utilities. These utilities should have no dependencies.*/
namespace chi_misc_utils
{
  std::string PrintIterationProgress(size_t current_iteration,
                                     size_t total_num_iterations,
                                     unsigned int num_intvls = 10);

//#include <iostream>
//#include <string>

const std::string WHITESPACE = " \n\r\t\f\v";

inline std::string ltrim(const std::string &s)
{
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

inline std::string rtrim(const std::string &s)
{
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

inline std::string trim(const std::string &s) {
  return rtrim(ltrim(s));
}
}//namespace chi_misc_utils

#endif //CHITECH_CHI_MISC_UTILS_H
