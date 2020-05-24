#ifndef _chi_volumemesher_predefined3d_h
#define _chi_volumemesher_predefined3d_h

#include "../chi_volumemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

class chi_mesh::VolumeMesherPredefined3D : public chi_mesh::VolumeMesher
{
public:
  void Execute();

  int GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid);
  bool IsRawCellNeighborToPartitionKBA(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell);
  void KBA(chi_mesh::UnpartitionedMesh* umesh,
           chi_mesh::MeshContinuum* grid);

  bool IsRawCellNeighborToPartitionParmetis(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    const std::vector<int>& cell_pids);
  void PARMETIS(chi_mesh::UnpartitionedMesh* umesh,
                chi_mesh::MeshContinuum* grid);
};


#endif