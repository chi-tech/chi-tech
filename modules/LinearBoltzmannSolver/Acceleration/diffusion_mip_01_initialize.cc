#include "diffusion_mip.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Initializes the diffusion solver. This involves creating the
 * sparse matrix with the appropriate sparsity pattern. Creating the
 * RHS vector. Creating the KSP solver. Setting the very specialized parameters
 * for Hypre's BooomerAMG. Note: `PCSetFromOptions` and
 * `KSPSetFromOptions` are called at the end. Therefore, any number of
 * additional PETSc options can be passed via the commandline.*/
void lbs::acceleration::DiffusionMIPSolver::Initialize()
{
  if (options.verbose)
    chi::log.Log() << m_text_name << ": Initializing PETSc items";

  if (options.verbose)
    chi::log.Log()
    << m_text_name
    << ": Global number of DOFs=" << m_num_global_dofs;

  //============================================= Create Matrix
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  m_sdm.BuildSparsityPattern(nodal_nnz_in_diag,
                             nodal_nnz_off_diag,
                             m_uk_man);

  m_A = chi_math::PETScUtils::CreateSquareMatrix(m_num_local_dofs,
                                                 m_num_global_dofs);
  chi_math::PETScUtils::InitMatrixSparsity(m_A,
                                           nodal_nnz_in_diag,
                                           nodal_nnz_off_diag);

  //============================================= Create RHS
  m_rhs = chi_math::PETScUtils::CreateVector(m_num_local_dofs,
                                             m_num_global_dofs);

  //============================================= Create KSP
  KSPCreate(PETSC_COMM_WORLD, &m_ksp);
  KSPSetOptionsPrefix(m_ksp, m_text_name.c_str());
  KSPSetType(m_ksp, KSPCG);

  KSPSetTolerances(m_ksp,1.e-50,
                   options.residual_tolerance,1.0e50,
                   options.max_iters);

  //============================================= Set Pre-conditioner
  PC pc;
  KSPGetPC(m_ksp,&pc);
//  PCSetType(pc, PCGAMG);
  PCSetType(pc,PCHYPRE);

  PCHYPRESetType(pc,"boomeramg");
  std::vector<std::string> pc_options =
    {"pc_hypre_boomeramg_agg_nl 1",
     "pc_hypre_boomeramg_P_max 4",
     "pc_hypre_boomeramg_grid_sweeps_coarse 1",
     "pc_hypre_boomeramg_max_levels 25",
     "pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi",
     "pc_hypre_boomeramg_coarsen_type HMIS",
     "pc_hypre_boomeramg_interp_type ext+i"};

  if (m_grid.Attributes() & chi_mesh::DIMENSION_2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  if (m_grid.Attributes() & chi_mesh::DIMENSION_3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.8");

  for (const auto& option : pc_options)
    PetscOptionsInsertString(nullptr, ("-"+m_text_name+option).c_str());

  PetscOptionsInsertString(nullptr, options.additional_options_string.c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(m_ksp);
}