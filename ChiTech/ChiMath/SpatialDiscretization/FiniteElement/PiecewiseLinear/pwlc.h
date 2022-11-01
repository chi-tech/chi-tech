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
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

//######################################################### Class def
namespace chi_math
{
  /**Generalization of the Galerkin Finite Element Method
     * with piecewise linear basis functions
     * for use by either a Continues Finite Element Method (CFEM)
     * or a Discontinuous Finite Element Method (DFEM). */
  class SpatialDiscretization_PWLC : public chi_math::SpatialDiscretization_FE
  {
  public:
    std::vector<std::shared_ptr<CellMappingFE_PWL>> cell_mappings;

  public:
    QuadratureLine          line_quad_order_arbitrary;
    QuadratureTriangle      tri_quad_order_arbitrary;
    QuadratureQuadrilateral quad_quad_order_arbitrary;
    QuadratureTetrahedron   tet_quad_order_arbitrary;
    QuadratureHexahedron    hex_quad_order_arbitrary;

    std::map<uint64_t, int64_t> node_mapping;

  //  std::vector<int> cell_local_block_address;
  //  std::vector<std::pair<int,int>> neighbor_cell_block_address;

  private:
  //  std::vector<chi_mesh::Cell*> neighbor_cells;
  //  std::vector<CellPWLFEValues*> neighbor_cell_fe_views;

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
        finite_element::SetupFlags setup_flags=
        finite_element::NO_FLAGS_SET,
        QuadratureOrder qorder =
        QuadratureOrder::SECOND,
        CoordinateSystemType in_cs_type =
        CoordinateSystemType::CARTESIAN)
    { if (in_grid == nullptr) throw std::invalid_argument(
        "Null supplied as grid to SpatialDiscretization_PWLC.");
      return std::shared_ptr<SpatialDiscretization_PWLC>(
      new SpatialDiscretization_PWLC(in_grid, setup_flags, qorder, in_cs_type));}

    //01
  private:
    std::shared_ptr<CellMappingFE_PWL> MakeCellMappingFE(const chi_mesh::Cell& cell) const;

  public:

    void PreComputeCellSDValues();
  //  void PreComputeNeighborCellSDValues(chi_mesh::MeshContinuumPtr grid);
    std::shared_ptr<CellMappingFE_PWL> GetCellMappingFE(uint64_t cell_local_index);

    void CreateCellMappings();

  private:
    //02
    void OrderNodes();

  public:
    //03
    void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                              std::vector<int64_t>& nodal_nnz_off_diag,
                              UnknownManager& unknown_manager) override;

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
    size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) override;
    size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) override;
  //  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
  //                               chi_math::UnknownManager* unknown_manager);
  //
  //  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
  //                                      chi_math::UnknownManager* unknown_manager,
  //                                      unsigned int unknown_id=0);

    size_t GetCellNumNodes(const chi_mesh::Cell& cell) const override
    {return cell.vertex_ids.size();}

    std::vector<chi_mesh::Vector3>
    GetCellNodeLocations(const chi_mesh::Cell& cell) const override
    {
      std::vector<chi_mesh::Vector3> node_locations;
      node_locations.reserve(cell.vertex_ids.size());

      for (auto& vid : cell.vertex_ids)
        node_locations.emplace_back(ref_grid->vertices[vid]);

      return node_locations;
    }

    void LocalizePETScVector(Vec petsc_vector,
                             std::vector<double>& local_vector,
                             UnknownManager& unknown_manager)
                             override;

    //FE-utils
    const finite_element::UnitIntegralData&
    GetUnitIntegrals(const chi_mesh::Cell& cell) override
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

    const finite_element::InternalQuadraturePointData&
    GetQPData_Volumetric(const chi_mesh::Cell& cell) override
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

    const finite_element::FaceQuadraturePointData&
    GetQPData_Surface(const chi_mesh::Cell& cell,
                      const unsigned int face) override
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
  };
}

#endif //SPATIAL_DISCRETIZATION_PWLC_H