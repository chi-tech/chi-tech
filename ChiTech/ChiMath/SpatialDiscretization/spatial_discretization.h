#ifndef SPATIAL_DISCRETIZATION_H
#define SPATIAL_DISCRETIZATION_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMath/Quadratures/quadrature.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMath/SpatialDiscretization/CellMappings/cell_mapping_base.h"

#include <petscksp.h>

#include <vector>
#include <map>

namespace chi_math
{
  class SpatialDiscretization
  {
  public:
    const SpatialDiscretizationType type;

    const chi_mesh::MeshContinuumConstPtr ref_grid;
    const CoordinateSystemType cs_type;

    const UnknownManager UNITARY_UNKNOWN_MANAGER;

  protected:
    std::vector<std::unique_ptr<CellMapping>> cell_mappings;
    std::map<uint64_t, std::shared_ptr<CellMapping>> nb_cell_mappings;

    uint64_t              local_block_address = 0;
    std::vector<uint64_t> locJ_block_address;
    std::vector<uint64_t> locJ_block_size;

    uint64_t local_base_block_size=0;
    uint64_t globl_base_block_size=0;

  protected:
    typedef SpatialDiscretizationType SDMType;
    //00
    explicit
    SpatialDiscretization(chi_mesh::MeshContinuumPtr& in_grid,
                          CoordinateSystemType in_cs_type,
                          SDMType in_type = SDMType::UNDEFINED) :
      type(in_type),
      ref_grid(in_grid),
      cs_type(in_cs_type),
      UNITARY_UNKNOWN_MANAGER({std::make_pair(chi_math::UnknownType::SCALAR,0)})
    {
    }

    //01 AddViewOfContinuum
  public:
    const CellMapping& GetCellMapping(const chi_mesh::Cell& cell) const;

    //02 OrderNodes

  public:
    //03
    virtual
    void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                              std::vector<int64_t>& nodal_nnz_off_diag,
                              const UnknownManager& unknown_manager) const = 0;

    //04 Mappings
    virtual
    int64_t MapDOF(const chi_mesh::Cell& cell,
                   unsigned int node,
                   const UnknownManager& unknown_manager,
                   unsigned int unknown_id,
                   unsigned int component) const =0;

    virtual
    int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                        unsigned int node,
                        const UnknownManager& unknown_manager,
                        unsigned int unknown_id,
                        unsigned int component) const =0;

    virtual
    int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const =0;
    virtual
    int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const =0;

    //05 Utils
    virtual
    size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const = 0;
    virtual
    size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const = 0;

    virtual
    size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const
    {
      return 0; //TODO: Remove
    }
    virtual
    std::vector<int64_t>
    GetGhostDOFIndices(const UnknownManager& unknown_manager) const
    {
      return {}; //TODO: Remove
    }

    virtual
    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const = 0;

    virtual
    std::vector<chi_mesh::Vector3>
      GetCellNodeLocations(const chi_mesh::Cell& cell) const = 0;

  public:
    /**Develops a localized view of a petsc vector.
     * Each spatial discretization has a specialization of this
     * method.*/
    virtual void LocalizePETScVector(Vec petsc_vector,
                                     std::vector<double>& local_vector,
                                     const UnknownManager& unknown_manager) const = 0;

  };
}

#endif