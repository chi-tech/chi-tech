#ifndef SPATIAL_DISCRETIZATION_PWLC_H
#define SPATIAL_DISCRETIZATION_PWLC_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/quadrature_hexahedron.h"

//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWLC : public SpatialDiscretization_FE
{
private:
  //00
  explicit SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr& in_grid,
      chi_math::finite_element::SetupFlags setup_flags,
      chi_math::QuadratureOrder qorder,
      chi_math::CoordinateSystemType in_cs_type);

  //01
  std::shared_ptr<CellMappingFE_PWL> MakeCellMappingFE(const chi_mesh::Cell& cell) const;

  void PreComputeCellSDValues() override;

  //02
  void OrderNodes();

public:
  unsigned int local_base_block_size=0;
  unsigned int globl_base_block_size=0;
  int local_block_address = 0;

  std::vector<std::shared_ptr<CellMappingFE_PWL>> cell_mappings;
  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;
  std::map<int,int> node_mapping;

  chi_math::QuadratureLine          line_quad_order_arbitrary;
  chi_math::QuadratureTriangle      tri_quad_order_arbitrary;
  chi_math::QuadratureQuadrilateral quad_quad_order_arbitrary;
  chi_math::QuadratureTetrahedron   tet_quad_order_arbitrary;
  chi_math::QuadratureHexahedron    hex_quad_order_arbitrary;

  //prevent anything else other than a shared pointer
  static std::shared_ptr<SpatialDiscretization_PWLC> New(
      chi_mesh::MeshContinuumPtr& in_grid,
      chi_math::finite_element::SetupFlags setup_flags=
      chi_math::finite_element::SetupFlags::NO_FLAGS_SET,
      chi_math::QuadratureOrder qorder =
      chi_math::QuadratureOrder::SECOND,
      chi_math::CoordinateSystemType in_cs_type =
      chi_math::CoordinateSystemType::CARTESIAN)
  { 
    if (in_grid == nullptr) throw std::invalid_argument(
        "Null supplied as grid to SpatialDiscretization_PWLC.");
    return std::shared_ptr<SpatialDiscretization_PWLC>(
    new SpatialDiscretization_PWLC(in_grid, setup_flags, qorder, in_cs_type));
  }

  //01
  std::shared_ptr<CellMappingFE_PWL> GetCellMappingFE(uint64_t cell_local_index);

  //03
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
      std::vector<int64_t>& nodal_nnz_off_diag,
      chi_math::UnknownManager& unknown_manager) override;

  //04 Mappings
  int64_t MapDOF(const chi_mesh::Cell& cell,
      unsigned int node,
      const chi_math::UnknownManager& unknown_manager,
      unsigned int unknown_id,
      unsigned int component) const override;

  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
      unsigned int node,
      const chi_math::UnknownManager& unknown_manager,
      unsigned int unknown_id,
      unsigned int component) const override;

  int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const override
  { 
    return MapDOF(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0);
  }

  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const override
  {
     return MapDOFLocal(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0);
  }

  //05
  size_t GetNumLocalDOFs(chi_math::UnknownManager& unknown_manager) override;

  size_t GetNumGlobalDOFs(chi_math::UnknownManager& unknown_manager) override;

  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
  {
    return cell.vertex_ids.size();
  }

  void LocalizePETScVector(Vec petsc_vector,
      std::vector<double>& local_vector,
      chi_math::UnknownManager& unknown_manager) override;

  //FE-utils
  const chi_math::finite_element::UnitIntegralData&
      GetUnitIntegrals(const chi_mesh::Cell& cell) override
  {
    return fe_unit_integrals.at(cell.local_id);
  }

  const chi_math::finite_element::InternalQuadraturePointData&
      GetQPData_Volumetric(const chi_mesh::Cell& cell) override
  {
    return fe_vol_qp_data.at(cell.local_id);
  }

  const chi_math::finite_element::FaceQuadraturePointData&
      GetQPData_Surface(const chi_mesh::Cell& cell, const unsigned int face) override
  {
    const auto& face_data = fe_srf_qp_data.at(cell.local_id);
    return face_data.at(face);
  }
};

#endif //SPATIAL_DISCRETIZATION_PWLC_H
