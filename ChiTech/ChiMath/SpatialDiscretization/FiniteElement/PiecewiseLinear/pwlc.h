#ifndef SPATIAL_DISCRETIZATION_PWLC_H
#define SPATIAL_DISCRETIZATION_PWLC_H

//#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "pwl_base.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/quadrature_hexahedron.h"

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
    std::map<uint64_t, int64_t> node_mapping;
    std::map<uint64_t, int64_t> m_ghost_node_mapping;

  //  std::vector<int> cell_local_block_address;
  //  std::vector<std::pair<int,int>> neighbor_cell_block_address;

  private:
    //00
    explicit
    SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr& in_grid,
                               finite_element::SetupFlags setup_flags,
                               QuadratureOrder qorder,
                               CoordinateSystemType in_cs_type);

  public:
    //prevent anything else other than a shared pointer
    static
    std::shared_ptr<SpatialDiscretization_PWLC>
    New(chi_mesh::MeshContinuumPtr& in_grid,
        finite_element::SetupFlags setup_flags = finite_element::NO_FLAGS_SET,
        QuadratureOrder qorder = QuadratureOrder::SECOND,
        CoordinateSystemType in_cs_type = CoordinateSystemType::CARTESIAN);

//    //01
//  public:
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

#endif //SPATIAL_DISCRETIZATION_PWLC_H