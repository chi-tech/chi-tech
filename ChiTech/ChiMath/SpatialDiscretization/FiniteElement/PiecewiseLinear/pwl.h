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

//    std::map<uint64_t, int64_t> node_mapping;

    std::vector<int64_t> cell_local_block_address;
    std::vector<std::pair<uint64_t, int64_t>> neighbor_cell_block_address;

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
        finite_element::SetupFlags setup_flags=
        finite_element::NO_FLAGS_SET,
        QuadratureOrder qorder =
        QuadratureOrder::SECOND,
        CoordinateSystemType in_cs_type =
        CoordinateSystemType::CARTESIAN);

//    //01
//  protected:
//    void PreComputeCellSDValues();
//    void PreComputeNeighborCellSDValues();
//
//    void CreateCellMappings();

  private:
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
    size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const override;
    size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const override;
    size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;
    std::vector<int64_t>
    GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override;

    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell) const override;

    void LocalizePETScVector(Vec petsc_vector,
                             std::vector<double>& local_vector,
                             const UnknownManager& unknown_manager)
                             const override;

    //FE-utils
    const finite_element::UnitIntegralData&
      GetUnitIntegrals(const chi_mesh::Cell& cell) override;

    const finite_element::InternalQuadraturePointData&
      GetQPData_Volumetric(const chi_mesh::Cell& cell) override;

    const finite_element::FaceQuadraturePointData&
      GetQPData_Surface(const chi_mesh::Cell& cell,
                        unsigned int face_index) override;
  };
}

#endif //SPATIAL_DISCRETIZATION_PWLD_H