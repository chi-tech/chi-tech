#ifndef SPATIAL_DISCRETIZATION_PWLD_H
#define SPATIAL_DISCRETIZATION_PWLD_H

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
class SpatialDiscretization_PWLD : public SpatialDiscretization_FE
{
public:
  std::vector<std::shared_ptr<CellMappingFE_PWL>> cell_mappings;

private:
  bool                     mapping_initialized=false;
  bool                     nb_mapping_initialized=false;
public:
  chi_math::QuadratureLine          line_quad_order_second;
  chi_math::QuadratureTriangle      tri_quad_order_second;
  chi_math::QuadratureQuadrilateral quad_quad_order_second;
  chi_math::QuadratureTetrahedron   tet_quad_order_second;
  chi_math::QuadratureHexahedron    hex_quad_order_second;

  chi_math::QuadratureLine          line_quad_order_arbitrary;
  chi_math::QuadratureTriangle      tri_quad_order_arbitrary;
  chi_math::QuadratureQuadrilateral quad_quad_order_arbitrary;
  chi_math::QuadratureTetrahedron   tet_quad_order_arbitrary;
  chi_math::QuadratureHexahedron    hex_quad_order_arbitrary;

//  std::map<int,int> node_mapping;

  int              local_block_address = 0;
  std::vector<int> cell_local_block_address;
  std::vector<std::pair<int,int>> neighbor_cell_block_address;

//  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

  unsigned int local_base_block_size=0;
  unsigned int globl_base_block_size=0;

private:
  std::map<uint64_t, chi_mesh::Cell*>  neighbor_cells;
  std::map<uint64_t, std::shared_ptr<CellMappingFE_PWL>> neighbor_cell_fe_views;

private:
  typedef chi_math::finite_element::UnitIntegralData UIData;
  typedef chi_math::finite_element::InternalQuadraturePointData QPDataVol;
  typedef chi_math::finite_element::FaceQuadraturePointData QPDataFace;

  std::map<uint64_t, UIData>                  nb_fe_unit_integrals;
  std::map<uint64_t, QPDataVol>               nb_fe_vol_qp_data;
  std::map<uint64_t, std::vector<QPDataFace>> nb_fe_srf_qp_data;

  bool nb_integral_data_initialized=false;
  bool nb_qp_data_initialized=false;

private:
  chi_math::finite_element::UnitIntegralData            scratch_intgl_data;
  chi_math::finite_element::InternalQuadraturePointData scratch_vol_qp_data;
  chi_math::finite_element::FaceQuadraturePointData     scratch_face_qp_data;

private:
  //00
  explicit
  SpatialDiscretization_PWLD(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder,
                             chi_math::CoordinateSystemType in_cs_type);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_PWLD>
  New(chi_mesh::MeshContinuumPtr& in_grid,
      chi_math::finite_element::SetupFlags setup_flags=
      chi_math::finite_element::SetupFlags::NO_FLAGS_SET,
      chi_math::QuadratureOrder qorder =
      chi_math::QuadratureOrder::SECOND,
      chi_math::CoordinateSystemType in_cs_type =
      chi_math::CoordinateSystemType::CARTESIAN)
  { if (in_grid == nullptr) throw std::invalid_argument(
      "Null supplied as grid to SpatialDiscretization_PWLD.");
    return std::shared_ptr<SpatialDiscretization_PWLD>(
    new SpatialDiscretization_PWLD(in_grid, setup_flags, qorder, in_cs_type));}

  //01
private:
  std::shared_ptr<CellMappingFE_PWL> MakeCellMappingFE(const chi_mesh::Cell& cell) const;

public:

  void PreComputeCellSDValues() override;
  void PreComputeNeighborCellSDValues();
  std::shared_ptr<CellMappingFE_PWL> GetCellMappingFE(int cell_local_index);
  chi_mesh::Cell&  GetNeighborCell(int cell_glob_index);
  std::shared_ptr<CellMappingFE_PWL> GetNeighborCellMappingFE(int cell_glob_index);

private:
  //02
  void OrderNodes();

public:
  //03
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;

  //04
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
  { return MapDOF(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }
  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const override
  { return MapDOFLocal(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }

  //05
  size_t GetNumLocalDOFs(chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_math::UnknownManager& unknown_manager) override;
//  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
//                               chi_math::UnknownManager* unknown_manager);
//
//  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
//                                      chi_math::UnknownManager* unknown_manager,
//                                      unsigned int unknown_id=0);

  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager& unknown_manager)
                           override;

  //FE-utils
  const chi_math::finite_element::UnitIntegralData&
    GetUnitIntegrals(const chi_mesh::Cell& cell) override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (integral_data_initialized)
        return fe_unit_integrals.at(cell.local_id);
      else
      {
        auto cell_fe_view = GetCellMappingFE(cell.local_id);
        scratch_intgl_data.Reset();
        cell_fe_view->ComputeUnitIntegrals(scratch_intgl_data);
        return scratch_intgl_data;
      }
    }
    else
    {
      if (nb_integral_data_initialized)
        return nb_fe_unit_integrals.at(cell.global_id);
      else
      {
        auto cell_fe_view = GetNeighborCellMappingFE(cell.global_id);
        cell_fe_view->ComputeUnitIntegrals(scratch_intgl_data);
        return scratch_intgl_data;
      }
    }
  }

  const chi_math::finite_element::InternalQuadraturePointData&
    GetQPData_Volumetric(const chi_mesh::Cell& cell) override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (qp_data_initialized)
        return fe_vol_qp_data.at(cell.local_id);
      else
      {
        auto cell_fe_view = GetCellMappingFE(cell.local_id);
        cell_fe_view->InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
        return scratch_vol_qp_data;
      }
    }
    else
    {
      if (nb_qp_data_initialized)
        return nb_fe_vol_qp_data.at(cell.global_id);
      else
      {
        auto cell_fe_view = GetNeighborCellMappingFE(cell.global_id);
        cell_fe_view->InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
        return scratch_vol_qp_data;
      }
    }
  }

  const chi_math::finite_element::FaceQuadraturePointData&
    GetQPData_Surface(const chi_mesh::Cell& cell,
                      const unsigned int face) override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (qp_data_initialized)
      {
        const auto& face_data = fe_srf_qp_data.at(cell.local_id);

        return face_data.at(face);
      }
      else
      {
        auto cell_fe_view = GetCellMappingFE(cell.local_id);
        cell_fe_view->InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
        return scratch_face_qp_data;
      }
    }
    else
    {
      if (nb_qp_data_initialized)
      {
        const auto& face_data = nb_fe_srf_qp_data.at(cell.global_id);

        return face_data.at(face);
      }
      else
      {
        auto cell_fe_view = GetNeighborCellMappingFE(cell.global_id);
        cell_fe_view->InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
        return scratch_face_qp_data;
      }
    }
  }
};

#endif //SPATIAL_DISCRETIZATION_PWLD_H