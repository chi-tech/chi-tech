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
    std::vector<CellFVValues*> cell_fv_views;

  private:
    std::vector<bool>  cell_view_added_flags;

  private:
    std::map<uint64_t, std::unique_ptr<chi_mesh::Cell>> neighbor_cells;
    std::map<uint64_t, CellFVValues*> neighbor_cell_fv_views;

  private:
    explicit
    SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr& in_grid,
                             CoordinateSystemType in_cs_type);

  public:
    virtual ~SpatialDiscretization_FV() = default;
    //prevent anything else other than a shared pointer
    static
    std::shared_ptr<SpatialDiscretization_FV>
    New(chi_mesh::MeshContinuumPtr& in_grid,
        CoordinateSystemType in_cs_type =
        CoordinateSystemType::CARTESIAN)
    { return std::shared_ptr<SpatialDiscretization_FV>(
      new SpatialDiscretization_FV(in_grid, in_cs_type));}

    //01
    void PreComputeCellSDValues();
    void PreComputeNeighborCellSDValues();

    CellFVValues* MapFeView(uint64_t cell_local_index);
    CellFVValues* MapNeighborFeView(uint64_t cell_global_index);

    //02 node ordering
  private:
    void OrderNodes();

    //03 sparsity
  public:
    void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                              std::vector<int64_t>& nodal_nnz_off_diag,
                              UnknownManager& unknown_manager) override;

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
    size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) override;
    size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) override;
    size_t GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
                           UnknownManager& unknown_manager);
    std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
                                        UnknownManager& unknown_manager,
                                        unsigned int unknown_id=0);

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
    {return 1;}

    std::vector<chi_mesh::Vector3>
      GetCellNodeLocations(const chi_mesh::Cell& cell) const override
    {
      std::vector<chi_mesh::Vector3> node_locations(1,cell.centroid);

      return node_locations;
    }

    void LocalizePETScVector(Vec petsc_vector,
                             std::vector<double>& local_vector,
                             UnknownManager& unknown_manager)
                             override;
  };
}


#endif