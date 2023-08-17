#ifndef CHITECH_MESHGENERATOR_H
#define CHITECH_MESHGENERATOR_H

#include "ChiObject.h"
#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

namespace chi
{
class GraphPartitioner;
}

namespace chi_mesh
{

class MeshContinuum;

/** Mesh generation can be very complicated in parallel. Some mesh formats
 * do not store connectivity information and therefore we have to establish
 * connectivity after the file is read. Most mesh formats are also not
 * partitioned for multiple processors when read/generated. The design of these
 * types of objects need to consider the previous generation of how we did this.
 * We had the concept of a VolumeMesher, which essentially took an
 * `UnpartitionedMesh` and converted it into a partitioned mesh (i.e. the
 * required `MeshContinuum` object). Up to now we've really only had two
 * variants. The `VolumeMesherPredefinedUnpartitioned` and the
 * `VolumeMesherExtruder`. Both of which operated on an `UnpartitionedMesh`
 * object. With this new design we want to unify these concepts to make them
 * more extendable and therefore we split a `MeshGenerator`'s execution into a
 * phase that generates an unpartitioned mesh and a phase that then converts
 * this mesh into real mesh (with both steps customizable). The phase that
 * creates the real mesh can be hooked up to a partitioner that can also be
 * designed to be pluggable.
 * */
class MeshGenerator : public ChiObject
{
public:
  /**Final execution step. */
  virtual void Execute();

  static chi::InputParameters GetInputParameters();
  explicit MeshGenerator(const chi::InputParameters& params);

protected:
  /**Virtual method to generate the unpartitioned mesh for the next step.*/
  virtual std::unique_ptr<UnpartitionedMesh>
  GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh);

  // 01
  /**Executes the partitioner and configures the mesh as a real mesh.*/
  virtual std::shared_ptr<MeshContinuum>
  SetupMesh(std::unique_ptr<UnpartitionedMesh> input_umesh_ptr);

  // 02 utils
  /**Determines if a cells needs to be included as a ghost or as a local cell.*/
  bool
  CellHasLocalScope(const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
                    uint64_t cell_global_id,
                    const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                    const std::vector<int64_t>& cell_partition_ids);

  /**Converts a light-weight cell to a real cell.*/
  static std::unique_ptr<chi_mesh::Cell>
  SetupCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
            uint64_t global_id,
            uint64_t partition_id,
            const std::vector<chi_mesh::Vector3>& vertices);

  const double scale_;
  const bool replicated_;
  std::vector<MeshGenerator*> inputs_;
  chi::GraphPartitioner* partitioner_ = nullptr;
};

} // namespace chi_mesh

#endif // CHITECH_MESHGENERATOR_H
