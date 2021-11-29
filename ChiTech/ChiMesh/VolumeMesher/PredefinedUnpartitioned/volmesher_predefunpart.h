#ifndef VOLMESHER_PREDEFUNPART_H
#define VOLMESHER_PREDEFUNPART_H

#include "../chi_volumemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

//###################################################################
/**This volume mesher merely applies a partitioning of an
 * unpartitioned mesh.*/
class chi_mesh::VolumeMesherPredefinedUnpartitioned :
                            public chi_mesh::VolumeMesher
{
public:
  void Execute() override;

  static
  bool CellHasLocalScope(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    uint64_t cell_global_id,
    const std::vector<std::set<uint64_t>>& vertex_subscriptions,
    const std::vector<int64_t>& cell_partition_ids);

  static
  std::vector<int64_t> KBA(const chi_mesh::UnpartitionedMesh& umesh);

  static
  std::vector<int64_t> PARMETIS(const UnpartitionedMesh &umesh);

  static chi_mesh::Cell* MakeCell(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    uint64_t global_id,
    uint64_t partition_id,
    const std::vector<chi_mesh::Vector3>& vertices);
};
#endif //VOLMESHER_PREDEFUNPART_H