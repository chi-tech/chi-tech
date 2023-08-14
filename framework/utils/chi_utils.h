#ifndef CHITECH_CHI_UTILS_H
#define CHITECH_CHI_UTILS_H

#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>


/**Miscellaneous utilities. These utilities should have no dependencies.*/
namespace chi
{
  std::string PrintIterationProgress(size_t current_iteration,
                                     size_t total_num_iterations,
                                     unsigned int num_intvls = 10);

//#include <iostream>
//#include <string>

  const std::string WHITESPACE = " \n\r\t\f\v";

  /**Trims whitespace from the front of a string.*/
  std::string StringLTrim(const std::string &s);

  /**Trims whitespace from the back of a string.*/
  std::string StringRTrim(const std::string &s);

  /**Trims whitespace from the front and back of a string.*/
  std::string StringTrim(const std::string &s);

  /**Splits a string using the given delimiter. Consecutive delimiters
   * are treated as one.*/
  std::vector<std::string> StringSplit(const std::string& input,
                                       const std::string& delim=" ");

  /**The string portion, from the rear of the input string, up to
   * encountering the search_string.*/
  std::string StringUpToFirstReverse(const std::string& input,
                                     const std::string& search_string);

  void AssertReadibleFile(const std::string& file_name);

  template<typename T, typename B>
  bool VectorListHas(const std::vector<T>& list, const B& val)
  {
    return std::find(list.begin(), list.end(), val) != list.end();
  }

  struct SubSetInfo
  {
    size_t ss_begin;
    size_t ss_end;
    size_t ss_size;
  };
  /**Subdivides a number of items (X) into a desired number of sub sets (Y). The
  * remainder of X/Y, i.e. r=X/Y obeys the indentity r < Y. These items will be
  * distributed to the first Y sub-sets. Example:
  * MakeSubSets(6659, 8) generates subsets of sizes
  * {833,833,833,832,832,832,832,832}.*/
  std::vector<SubSetInfo> MakeSubSets(size_t num_items,
                                      size_t desired_num_subsets);
}//namespace chi_misc_utils

#endif // CHITECH_CHI_UTILS_H
