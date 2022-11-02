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
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

//######################################################### Class def
namespace chi_math
{
  /**Generalization of the Galerkin Finite Element Method
     * with piecewise linear basis functions
     * for use by either a Continues Finite Element Method (CFEM)
     * or a Discontinuous Finite Element Method (DFEM). */
  class SpatialDiscretization_PWLD : public chi_math::SpatialDiscretization_FE
  {
  protected:
    QuadratureLine          line_quad_order_arbitrary;
    QuadratureTriangle      tri_quad_order_arbitrary;
    QuadratureQuadrilateral quad_quad_order_arbitrary;
    QuadratureTetrahedron   tet_quad_order_arbitrary;
    QuadratureHexahedron    hex_quad_order_arbitrary;

    std::map<uint64_t, int64_t> node_mapping;

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
        CoordinateSystemType::CARTESIAN)
    { if (in_grid == nullptr) throw std::invalid_argument(
        "Null supplied as grid to SpatialDiscretization_PWLD.");
      return std::shared_ptr<SpatialDiscretization_PWLD>(
      new SpatialDiscretization_PWLD(in_grid, setup_flags, qorder, in_cs_type));}

    //01
  protected:
    void PreComputeCellSDValues();
    void PreComputeNeighborCellSDValues();

    void CreateCellMappings();

  private:
    //02
    void OrderNodes();

  public:
    //03
    void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                              std::vector<int64_t>& nodal_nnz_off_diag,
                              UnknownManager& unknown_manager) override;

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
      if (ref_grid->IsCellLocal(cell.global_id))
      {
        if (integral_data_initialized)
          return fe_unit_integrals.at(cell.local_id);
        else
        {
          const auto& cell_mapping = GetCellMapping(cell);
          scratch_intgl_data.Reset();
          cell_mapping.ComputeUnitIntegrals(scratch_intgl_data);
          return scratch_intgl_data;
        }
      }
      else
      {
        if (nb_integral_data_initialized)
          return nb_fe_unit_integrals.at(cell.global_id);
        else
        {
          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.ComputeUnitIntegrals(scratch_intgl_data);
          return scratch_intgl_data;
        }
      }
    }

    const finite_element::InternalQuadraturePointData&
      GetQPData_Volumetric(const chi_mesh::Cell& cell) override
    {
      if (ref_grid->IsCellLocal(cell.global_id))
      {
        if (qp_data_initialized)
          return fe_vol_qp_data.at(cell.local_id);
        else
        {
          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
          return scratch_vol_qp_data;
        }
      }
      else
      {
        if (nb_qp_data_initialized)
          return nb_fe_vol_qp_data.at(cell.global_id);
        else
        {
          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
          return scratch_vol_qp_data;
        }
      }
    }

    const finite_element::FaceQuadraturePointData&
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
          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
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
          const auto& cell_mapping = GetCellMapping(cell);
          cell_mapping.InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
          return scratch_face_qp_data;
        }
      }
    }
  };
}

#endif //SPATIAL_DISCRETIZATION_PWLD_H