#ifndef SPATIAL_DISCRETIZATION_FV_H
#define SPATIAL_DISCRETIZATION_FV_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "CellViews/fv_cellbase.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "CellViews/fv_cellbase.h"

#include <map>

//###################################################################
namespace chi_math
{
/**Spatial discretizations supporting Finite Volume representations.
   * */
class SpatialDiscretization_FV : public chi_math::SpatialDiscretization
{
private:
  std::map<uint64_t, uint64_t> neighbor_cell_local_ids_;
private:
  explicit
  SpatialDiscretization_FV(const chi_mesh::MeshContinuum& in_grid,
                           CoordinateSystemType in_cs_type);

public:
  virtual ~SpatialDiscretization_FV() = default;
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_FV>
  New(const chi_mesh::MeshContinuum& in_grid,
      CoordinateSystemType in_cs_type =
      CoordinateSystemType::CARTESIAN)
  { return std::shared_ptr<SpatialDiscretization_FV>(
    new SpatialDiscretization_FV(in_grid, in_cs_type));}

  //01
  void CreateCellMappings();

  //02 node ordering
protected:
  void OrderNodes();

  //03 sparsity
public:
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            const UnknownManager& unknown_manager) const override;

  //04a mappings
  int64_t MapDOF(const chi_mesh::Cell& cell,
                 unsigned int node,
                 const UnknownManager& unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component) const override;

  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                  unsigned int node,
                  const UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component) const override;

  int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const override
  { return MapDOF(cell,node,UNITARY_UNKNOWN_MANAGER,0,0); }
  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const override
  { return MapDOFLocal(cell,node,UNITARY_UNKNOWN_MANAGER,0,0); }

  //05 utils
  size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const override;
  size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const override;
  size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;
  std::vector<int64_t>
    GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;

  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
  {return 1;}

  std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell) const override
  {
    std::vector<chi_mesh::Vector3> node_locations(1,cell.centroid_);

    return node_locations;
  }
};

}//namespace chi_math


#endif