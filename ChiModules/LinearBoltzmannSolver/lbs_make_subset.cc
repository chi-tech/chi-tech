#include "lbs_make_subset.h"

#include <cmath>

namespace lbs
{
std::vector<SubSetInfo> MakeSubSets(const size_t num_items,
                                    const size_t desired_num_subsets)
{
  //=================================== Set number of subsets
  //We cannot have more subsets than there
  //are items
  const size_t num_subsets = std::min(desired_num_subsets, num_items);

  //=================================== Determine subset size
  //This is the general size of each subset
  //except the last. We use the floor function
  //here because we want the last set to contain
  //the remainder
  const size_t subset_gen_size = std::floor(num_items / num_subsets);

  //=================================== Populate the subsets
  std::vector<SubSetInfo> ss_infos;
  ss_infos.reserve(num_subsets);
  for (size_t ssid=0; ssid<num_subsets; ++ssid)
  {
    const size_t ss_start = ssid * subset_gen_size;
    const size_t ss_end   = (ssid < (num_subsets-1))?
                            (ssid+1) * subset_gen_size-1 :
                            num_items-1;
    const size_t ss_size  = ss_end-ss_start+1;

    ss_infos.push_back({ss_start, ss_end, ss_size});
  }

  return ss_infos;
}

}//namespace lbs

