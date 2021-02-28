#ifndef _volmesher_predefunpart_h
#define _volmesher_predefunpart_h

#include "../chi_volumemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

//###################################################################
/**This volume mesher merely applies a partitioning of an
 * unpartitioned mesh.*/
class chi_mesh::VolumeMesherPredefinedUnpartitioned :
                            public chi_mesh::VolumeMesher
{
public:
  void Execute();

  void BuildMeshConnectivity(chi_mesh::UnpartitionedMesh* umesh);

  int GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid);

  bool IsRawCellNeighborToPartitionKBA(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell);

  void KBA(chi_mesh::UnpartitionedMesh* umesh,
           chi_mesh::MeshContinuumPtr grid);

  bool IsRawCellNeighborToPartitionParmetis(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    const std::vector<int>& cell_pids);

  void PARMETIS(chi_mesh::UnpartitionedMesh* umesh,
                chi_mesh::MeshContinuumPtr grid);

  void AddSlabToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
  void AddPolygonToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
  void AddPolyhedronToGrid(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& raw_cell,
    const chi_mesh::Cell& temp_cell,
    chi_mesh::MeshContinuum& grid);
};
#endif