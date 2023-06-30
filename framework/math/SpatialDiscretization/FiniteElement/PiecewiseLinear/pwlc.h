#ifndef SPATIAL_DISCRETIZATION_PWLC_H
#define SPATIAL_DISCRETIZATION_PWLC_H

#include "pwl_base.h"
#include "math/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

//######################################################### Class def
namespace chi_math
{
  /**Generalization of the Galerkin Finite Element Method
     * with piecewise linear basis functions
     * for use by either a Continues Finite Element Method (CFEM)
     * or a Discontinuous Finite Element Method (DFEM). */
  class SpatialDiscretization_PWLC : public SpatialDiscretization_PWLBase
  {
  protected:
    std::map<uint64_t, int64_t> node_mapping_;
    std::map<uint64_t, int64_t> ghost_node_mapping_;

  private:
    //00
    explicit
    SpatialDiscretization_PWLC(const chi_mesh::MeshContinuum& in_grid,
                               finite_element::SetupFlags setup_flags,
                               QuadratureOrder qorder,
                               CoordinateSystemType in_cs_type);

  public:
    //prevent anything else other than a shared pointer
    static
    std::shared_ptr<SpatialDiscretization_PWLC>
    New(const chi_mesh::MeshContinuum& in_grid,
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

    //04 Mappings
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

#endif //SPATIAL_DISCRETIZATION_PWLC_H