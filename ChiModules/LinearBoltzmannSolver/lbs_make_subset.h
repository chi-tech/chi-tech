#ifndef CHITECH_LBS_MAKE_SUBSET_H
#define CHITECH_LBS_MAKE_SUBSET_H

#include <vector>
#include <cstddef>

namespace lbs
{

struct SubSetInfo
{
  size_t ss_begin;
  size_t ss_end;
  size_t ss_size;
};

std::vector<SubSetInfo> MakeSubSets(size_t num_items,
                                    size_t desired_num_subsets);


}//namespace lbs

#endif //CHITECH_LBS_MAKE_SUBSET_H
