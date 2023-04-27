#include "diffusion.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

namespace lbs::acceleration
{

// ###################################################################
/**Default constructor.*/
DiffusionSolver::DiffusionSolver(
  std::string text_name,
  const chi_math::SpatialDiscretization& sdm,
  const chi_math::UnknownManager& uk_man,
  std::map<uint64_t, BoundaryCondition> bcs,
  MatID2XSMap map_mat_id_2_xs,
  const std::vector<UnitCellMatrices>& unit_cell_matrices,
  const bool verbose,
  const bool requires_ghosts)
  : text_name_(std::move(text_name)),
    grid_(sdm.ref_grid_),
    sdm_(sdm),
    uk_man_(uk_man),
    bcs_(std::move(bcs)),
    mat_id_2_xs_map_(std::move(map_mat_id_2_xs)),
    unit_cell_matrices_(unit_cell_matrices),
    num_local_dofs_(static_cast<int64_t>(sdm_.GetNumLocalDOFs(uk_man_))),
    num_global_dofs_(static_cast<int64_t>(sdm_.GetNumGlobalDOFs(uk_man_))),
    A_(nullptr),
    rhs_(nullptr),
    ksp_(nullptr),
    requires_ghosts_(requires_ghosts)
{
  options.verbose = verbose;
}

// ###################################################################
/**Default destructor.*/
DiffusionSolver::~DiffusionSolver()
{
  MatDestroy(&A_);
  VecDestroy(&rhs_);
  KSPDestroy(&ksp_);
}

// ###################################################################
/**Returns the assigned text name.*/
std::string DiffusionSolver::TextName() const {return text_name_;}

// ###################################################################
/**Returns the right-hand side petsc vector.*/
const Vec& DiffusionSolver::RHS() const { return rhs_; }

// ###################################################################
/**Returns the assigned unknown structure.*/
const chi_math::UnknownManager& DiffusionSolver::UnknownStructure() const
{
  return uk_man_;
}

// ###################################################################
/**Returns the associated spatial discretization.*/
const chi_math::SpatialDiscretization&
DiffusionSolver::SpatialDiscretization() const
{
  return sdm_;
}

std::pair<size_t, size_t>
lbs::acceleration::DiffusionSolver::GetNumPhiIterativeUnknowns()
{
  return {sdm_.GetNumLocalDOFs(uk_man_), sdm_.GetNumGlobalDOFs(uk_man_)};
}



} // namespace lbs::acceleration
