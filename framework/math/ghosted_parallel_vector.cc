#include "ghosted_parallel_vector.h"

#include "chi_log_exceptions.h"

#include <algorithm>


namespace chi_math
{

double GhostedParallelVector::GetGlobalValue(const int64_t global_id) const
{
  if (global_id >= extents_[location_id_] and
      global_id < extents_[location_id_ + 1])
    return values_[global_id - extents_[location_id_]];

  ChiInvalidArgumentIf(
      ghost_comm_.MapGhostToLocal(global_id) == -1,
      "Invalid global id specified. Specified global ids must be "
      "locally owned or ghosts.");
  return values_[ghost_comm_.MapGhostToLocal(global_id)];
}


}
