#ifndef CHI_DISCRETIZATION_PWL_H
#define CHI_DISCRETIZATION_PWL_H

#include"ChiMath/SpatialDiscretization/spatial_discretization.h"
#include"../../../ChiMesh/Region/chi_region.h"
#include "CellViews/pwl_cellbase.h"
#include "../../Quadratures/quadrature_triangle.h"
#include "../../Quadratures/quadrature_tetrahedron.h"





//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWL : public SpatialDiscretization
{
public:
  std::vector<CellFEValues*> cell_fe_views;

private:
  std::vector<bool>        cell_view_added_flags;
  bool                     mapping_initialized;
public:
  chi_math::QuadratureTriangle*    tri_quad_deg5;
  chi_math::QuadratureTriangle*    tri_quad_deg3_surf;
  chi_math::QuadratureTetrahedron* tet_quad_deg1;
  chi_math::QuadratureTetrahedron* tet_quad_deg3;
  chi_math::QuadratureTetrahedron* tet_quad_deg3_surface;

private:
  std::vector<chi_mesh::Cell*> neighbor_cells;
  std::vector<CellFEValues*> neighbor_cell_fe_views;

  typedef chi_math::SpatialDiscretizationType SDMType;

private:
  //00
  explicit
  SpatialDiscretization_PWL(int dim=0, SDMType sd_method =
                                       SDMType::PIECEWISE_LINEAR_DISCONTINUOUS);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_PWL>
  New(int in_dim=0, SDMType in_sd_method =
                    SDMType::PIECEWISE_LINEAR_DISCONTINUOUS)
  { return std::shared_ptr<SpatialDiscretization_PWL>(
    new SpatialDiscretization_PWL(in_dim, in_sd_method));}

  //01
  void PreComputeCellSDValues(chi_mesh::MeshContinuumPtr grid) override;
  void AddViewOfNeighborContinuums(chi_mesh::MeshContinuumPtr grid);
  CellFEValues* MapFeViewL(int cell_local_index);

  //02
  std::pair<int,int> OrderNodesCFEM(chi_mesh::MeshContinuumPtr grid);

  //03
  void BuildCFEMSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                const std::pair<int,int>& domain_ownership);
  void BuildCFEMSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                chi_math::NodalVariableStructure* unknown_manager=nullptr);

  //04
  std::pair<int,int> OrderNodesDFEM(chi_mesh::MeshContinuumPtr grid);

  //05
  void BuildDFEMSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                const std::pair<int,int>& domain_ownership);
  void BuildDFEMSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                                std::vector<int>& nodal_nnz_in_diag,
                                std::vector<int>& nodal_nnz_off_diag,
                                chi_math::NodalVariableStructure* unknown_manager=nullptr);
  chi_mesh::Cell* MapNeighborCell(int cell_glob_index);
  CellFEValues* MapNeighborCellFeView(int cell_glob_index);

  //06a Mappings
  int MapCFEMDOF(int vertex_id);
  int MapCFEMDOF(int vertex_id,
                 chi_math::NodalVariableStructure* unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component=0);

  int MapDFEMDOF(chi_mesh::Cell* cell, int dof,
                 int component=0,
                 int component_block_offset=1);
  int MapDFEMDOFLocal(chi_mesh::Cell* cell, int dof,
                      int component=0,
                      int component_block_offset=1);

  int MapDFEMDOF(chi_mesh::Cell* cell, int node,
                 chi_math::NodalVariableStructure* unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component=0);
  int MapDFEMDOFLocal(chi_mesh::Cell* cell, int node,
                      chi_math::NodalVariableStructure* unknown_manager,
                      unsigned int unknown_id,
                      unsigned int component=0);


  //06b utils
  unsigned int GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                               chi_math::NodalVariableStructure* unknown_manager=nullptr);
  unsigned int GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                                chi_math::NodalVariableStructure* unknown_manager=nullptr);
  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
                               chi_math::NodalVariableStructure* unknown_manager);

  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
                                      chi_math::NodalVariableStructure* unknown_manager,
                                      unsigned int unknown_id=0);

  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::NodalVariableStructure* unknown_manager)
                           override;
};

#endif