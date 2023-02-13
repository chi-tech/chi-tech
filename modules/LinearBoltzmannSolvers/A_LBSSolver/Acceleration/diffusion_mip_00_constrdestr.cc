#include <utility>

#include "acceleration.h"
#include "diffusion_mip.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

//###################################################################
/**Default constructor.*/
lbs::acceleration::DiffusionMIPSolver::
  DiffusionMIPSolver(std::string text_name,
                     const chi_mesh::MeshContinuum &grid,
                     const chi_math::SpatialDiscretization& sdm,
                     const chi_math::UnknownManager& uk_man,
                     std::map<uint64_t, BoundaryCondition>  bcs,
                     MapMatID2XS map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     const bool verbose/*=false*/) :
  m_text_name(std::move(text_name)),
  m_grid(grid),
  m_sdm(sdm),
  m_uk_man(uk_man),
  m_bcs(std::move(bcs)),
  m_map_mat_id_2_xs(std::move(map_mat_id_2_xs)),
  m_unit_cell_matrices(unit_cell_matrices),
  m_num_local_dofs(static_cast<int64_t>(m_sdm.GetNumLocalDOFs(m_uk_man))),
  m_num_global_dofs(static_cast<int64_t>(m_sdm.GetNumGlobalDOFs(m_uk_man))),
  m_A(nullptr),
  m_rhs(nullptr),
  m_ksp(nullptr)
{
  options.verbose = verbose;

  using SDM_TYPE = chi_math::SpatialDiscretizationType;
  const auto& PWLD = SDM_TYPE ::PIECEWISE_LINEAR_DISCONTINUOUS;

  if (m_sdm.type != PWLD)
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver: can only be"
                           " used with PWLD.");
}

//###################################################################
/**Default destructor.*/
lbs::acceleration::DiffusionMIPSolver::~DiffusionMIPSolver()
{
  MatDestroy(&m_A);
  VecDestroy(&m_rhs);
  KSPDestroy(&m_ksp);
}

const Vec& lbs::acceleration::DiffusionMIPSolver::RHS() const
{
  return m_rhs;
}