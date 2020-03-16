#ifndef _chi_volumemesher_predefined3d_h
#define _chi_volumemesher_predefined3d_h

#include "../chi_volumemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

class chi_mesh::VolumeMesherPredefined3D : public chi_mesh::VolumeMesher
{
public:
  int GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid);
  bool IsRawCellNeighborToPartition(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell);
  void Execute();
};


#endif