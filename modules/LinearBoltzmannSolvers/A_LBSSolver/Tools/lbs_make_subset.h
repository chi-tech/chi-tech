#ifndef CHITECH_LBS_MAKE_SUBSET_H
#define CHITECH_LBS_MAKE_SUBSET_H

#include <vector>
#include <cstddef>

namespace lbs
{

/**Data structure for information related to a subset of items.*/
struct SubSetInfo
{
  size_t ss_begin;
  size_t ss_end;
  size_t ss_size;
};

/**Routine that neatly contains the logic for creating a number of
 * subsets from a given set.*/
std::vector<SubSetInfo> MakeSubSets(size_t num_items,
                                    size_t desired_num_subsets);


}//namespace lbs

#endif //CHITECH_LBS_MAKE_SUBSET_H
