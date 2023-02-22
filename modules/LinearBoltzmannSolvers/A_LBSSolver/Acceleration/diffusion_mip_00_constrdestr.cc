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
                     MatID2XSMap map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     const bool verbose/*=false*/) :
  text_name_(std::move(text_name)),
  grid_(grid),
  sdm_(sdm),
  uk_man_(uk_man),
  bcs_(std::move(bcs)),
  mat_id_2_xs_map(std::move(map_mat_id_2_xs)),
  unit_cell_matrices_(unit_cell_matrices),
  num_local_dofs_(static_cast<int64_t>(sdm_.GetNumLocalDOFs(uk_man_))),
  num_global_dofs_(static_cast<int64_t>(sdm_.GetNumGlobalDOFs(uk_man_))),
  A_(nullptr),
  rhs_(nullptr),
  ksp_(nullptr)
{
  options.verbose = verbose;

  using SDM_TYPE = chi_math::SpatialDiscretizationType;
  const auto& PWLD = SDM_TYPE ::PIECEWISE_LINEAR_DISCONTINUOUS;

  if (sdm_.type != PWLD)
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver: can only be"
                           " used with PWLD.");
}

//###################################################################
/**Default destructor.*/
lbs::acceleration::DiffusionMIPSolver::~DiffusionMIPSolver()
{
  MatDestroy(&A_);
  VecDestroy(&rhs_);
  KSPDestroy(&ksp_);
}

const Vec& lbs::acceleration::DiffusionMIPSolver::RHS() const
{
  return rhs_;
}