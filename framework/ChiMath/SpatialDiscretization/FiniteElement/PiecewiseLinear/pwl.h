#ifndef SPATIAL_DISCRETIZATION_PWLD_H
#define SPATIAL_DISCRETIZATION_PWLD_H

#include "pwl_base.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

//######################################################### Class def
namespace chi_math
{
  /**Generalization of the Galerkin Finite Element Method
     * with piecewise linear basis functions
     * for use by either a Continues Finite Element Method (CFEM)
     * or a Discontinuous Finite Element Method (DFEM). */
  class SpatialDiscretization_PWLD : public SpatialDiscretization_PWLBase
  {
  protected:
    std::vector<int64_t> cell_local_block_address_;
    std::vector<std::pair<uint64_t, int64_t>> neighbor_cell_block_address_;

  private:
    //00
    explicit
    SpatialDiscretization_PWLD(chi_mesh::MeshContinuumPtr& in_grid,
                               finite_element::SetupFlags setup_flags,
                               QuadratureOrder qorder,
                               CoordinateSystemType in_cs_type);

  public:
    //prevent anything else other than a shared pointer
    static
    std::shared_ptr<SpatialDiscretization_PWLD>
    New(chi_mesh::MeshContinuumPtr& in_grid,
        finite_element::SetupFlags setup_flags = finite_element::NO_FLAGS_SET,
        QuadratureOrder qorder = QuadratureOrder::SECOND,
        CoordinateSystemType in_cs_type = CoordinateSystemType::CARTESIAN);

    //01
    //Inherited from PWLBase:
    //PreComputeCellSDValues
    //PreComputeNeighborCellSDValues
    //CreateCellMappings

  protected:
    //02
    void OrderNodes();

  public:
    //03
    void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                              std::vector<int64_t>& nodal_nnz_off_diag,
                              const UnknownManager& unknown_manager) const override;

    //04
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

    //05
    //Inherited from PWLBase:
    //GetNumLocalDOFs
    //GetNumGlobalDOFs

    size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;

    std::vector<int64_t>
    GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;

    //Inherited from PWLBase:
    //GetCellNumNodes
    //GetCellNodeLocations
    //LocalizePETScVector
    //GetUnitIntegrals
    //GetQPData_Volumetric
    //GetQPData_Surface
  };
}

#endif //SPATIAL_DISCRETIZATION_PWLD_H