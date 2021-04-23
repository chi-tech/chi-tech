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
  int GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid);

  static
  bool IsRawCellNeighborToPartitionKBA(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell);

  static
  void KBA(chi_mesh::UnpartitionedMesh* umesh,
           chi_mesh::MeshContinuumPtr& grid);

  static
  bool IsRawCellNeighborToPartitionParmetis(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    const std::vector<int64_t>& cell_pids);

  static
  void PARMETIS(chi_mesh::UnpartitionedMesh* umesh,
                chi_mesh::MeshContinuumPtr& grid);

  static
  void AddSlabToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
  static
  void AddPolygonToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
  static
  void AddPolyhedronToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
};
#endif //VOLMESHER_PREDEFUNPART_H