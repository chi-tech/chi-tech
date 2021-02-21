#ifndef SPATIAL_DISCRETIZATION_PWLD_H
#define SPATIAL_DISCRETIZATION_PWLD_H

#include "CHI_TECH/ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_cellbase.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"

//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWL : public SpatialDiscretization_FE
{
public:
  std::vector<CellPWLFEValues*> cell_fe_views;

private:
  bool                     mapping_initialized=false;
public:
  chi_math::QuadratureGaussLegendre line_quad_order_second;
  chi_math::QuadratureTriangle      tri_quad_order_second;
  chi_math::QuadratureTetrahedron   tet_quad_order_second;

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
  std::map<uint64_t, CellPWLFEValues*> neighbor_cell_fe_views;

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
  //00
  explicit
  SpatialDiscretization_PWL(chi_mesh::MeshContinuumPtr in_grid);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_PWL>
  New(chi_mesh::MeshContinuumPtr in_grid)
  { if (in_grid == nullptr) throw std::invalid_argument(
        "Null supplied as grid to SpatialDiscretization_PWLC.");
    return std::shared_ptr<SpatialDiscretization_PWL>(
    new SpatialDiscretization_PWL(in_grid));}

  //01
  void PreComputeCellSDValues(chi_mesh::MeshContinuumPtr grid) override;
  void PreComputeNeighborCellSDValues(chi_mesh::MeshContinuumPtr grid);
  CellPWLFEValues* MapFeViewL(int cell_local_index);

private:
  //02
  void OrderNodes(chi_mesh::MeshContinuumPtr grid);

public:
  //03
  void BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;
  chi_mesh::Cell* MapNeighborCell(int cell_glob_index);
  CellPWLFEValues* MapNeighborCellFeView(int cell_glob_index);

  //04
  int MapDOF(const chi_mesh::Cell& cell,
             const int node,
             const chi_math::UnknownManager& unknown_manager,
             const unsigned int unknown_id,
             const unsigned int component= 0) const;
  int MapDOFLocal(const chi_mesh::Cell& cell,
                  const int node,
                  const chi_math::UnknownManager& unknown_manager,
                  const unsigned int unknown_id,
                  const unsigned int component= 0) const;
  int MapDOF(const chi_mesh::Cell& cell, int node)
  { return MapDOF(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }
  int MapDOFLocal(const chi_mesh::Cell& cell, int node)
  { return MapDOFLocal(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }

  //05
  size_t GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                          chi_math::UnknownManager& unknown_manager) override;
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
    GetUnitIntegrals(const chi_mesh::Cell& cell) const override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (not integral_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetUnitIntegrals "
                                    "called without integrals being initialized.");
      return fe_unit_integrals[cell.local_id];
    }
    else
    {
      if (not nb_integral_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetUnitIntegrals "
                                    "called without integrals being initialized.");
      return nb_fe_unit_integrals.at(cell.global_id);
    }
  }

  const chi_math::finite_element::InternalQuadraturePointData&
    GetQPData_Volumetric(const chi_mesh::Cell& cell) const override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (not qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetQPData_Volumetric "
                                    "called without integrals being initialized.");
      return fe_vol_qp_data.at(cell.local_id);
    }
    else
    {
      if (not nb_qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetQPData_Volumetric "
                                    "called without quadrature data being initialized.");
      return nb_fe_vol_qp_data.at(cell.global_id);
    }
  }

  const chi_math::finite_element::FaceQuadraturePointData&
    GetQPData_Surface(const chi_mesh::Cell& cell,
                      const unsigned int face) const override
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      if (not qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetQPData_Surface "
                                    "called without quadrature data being initialized.");

      const auto& face_data = fe_srf_qp_data.at(cell.local_id);

      return face_data.at(face);
    }
    else
    {
      if (not nb_qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_PWLD::GetQPData_Volumetric "
                                    "called without quadrature data being initialized.");

      const auto& face_data = nb_fe_srf_qp_data.at(cell.global_id);

      return face_data.at(face);
    }
  }
};

#endif //SPATIAL_DISCRETIZATION_PWLD_H