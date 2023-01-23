#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//============================================= constructor
mg_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, { {"max_inner_iters"   , int64_t(500)},
                                        {"residual_tolerance", 1.0e-2},
                                        {"verbose_level"     , int64_t (0) },
                                        {"thermal_flux_tolerance", 1.0e-2},
                                        {"max_thermal_iters" , int64_t(500)},
                                        {"do_two_grid"       , false}
  })
{}

//============================================= destructor
mg_diffusion::Solver::~Solver()
{
  for (uint g=0; g<num_groups; ++g)
  {
    VecDestroy(&x[g]);
    VecDestroy(&bext[g]);
    MatDestroy(&A[g]);
  }
  VecDestroy(&b);

  if (last_fast_group < num_groups)
  {
    VecDestroy(&thermal_dphi);
    for (uint g=last_fast_group; g<num_groups; ++g)
      VecDestroy(&x_old[g]);
  }

  if (do_two_grid)
  {
    VecDestroy(&x[num_groups]);
    MatDestroy(&A[num_groups]);
  }
}

//============================================= Initialize
void mg_diffusion::Solver::Initialize()
{
  chi::log.Log() << "\n"
                 << chi::program_timer.GetTimeString() << " "
                 << TextName() << ": Initializing CFEM Multigroup Diffusion solver ";

  //============================================= Get grid
  grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;
  if (grid_ptr == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " No grid defined.");
 
  chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //================================================== Add unique material ids
  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (auto& cell : grid.local_cells)
  {
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id<0)
      ++invalid_mat_cell_count;
  }
  const auto& ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t cell_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[cell_id];
    unique_material_ids.insert(cell.material_id);
    if (cell.material_id<0)
      ++invalid_mat_cell_count;
  }

  if (invalid_mat_cell_count>0)
  {
    chi::log.LogAllWarning()
      << "Number of invalid material cells: " << invalid_mat_cell_count;
  }

  //================================================== Initialize materials
  mg_diffusion::Solver::Initialize_Materials(unique_material_ids);

  //============================================= BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();
  mg_diffusion::Solver::Set_BCs(globl_unique_bndry_ids);
  
  //============================================= Make SDM
  sdm_ptr = chi_math::SpatialDiscretization_PWLC::New(grid_ptr);
  const auto& sdm = *sdm_ptr;
 
  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);
 
  chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag, OneDofPerNode);

  unsigned int i_two_grid = do_two_grid ? 1 : 0;
  std::cout << "i_two_grid = " << i_two_grid << std::endl;

  A.resize(num_groups+i_two_grid, nullptr);
  x.resize(num_groups+i_two_grid, nullptr);
  bext.resize(num_groups, nullptr);

  auto ghost_dof_indices = sdm.GetGhostDOFIndices(OneDofPerNode);

  for (uint g=0; g < num_groups; ++g)
  {
    // x[g] = chi_math::PETScUtils::CreateVector(n,N);
    x[g] = chi_math::PETScUtils::CreateVectorWithGhosts(n,N,
                  static_cast<int64_t>(ghost_dof_indices.size()),
                  ghost_dof_indices);
    VecSet(x[g], 0.0);
    bext[g] = chi_math::PETScUtils::CreateVector(n,N);

    A[g] = chi_math::PETScUtils::CreateSquareMatrix(n,N);
    chi_math::PETScUtils::InitMatrixSparsity(A[g],
                                             nodal_nnz_in_diag,
                                             nodal_nnz_off_diag);
  }
  // initialize b
  VecDuplicate(bext.front(), &b);
  // initialize old flux for thermal groups
  if (last_fast_group < num_groups)
  {
    VecDuplicate(x.front(), &thermal_dphi);
    x_old.resize(num_groups, nullptr);
    for (uint g=0; g < num_groups; ++g)
    {
      VecDuplicate(x.front(), &x_old[g]); //jcr is x_old like bext or like x?
      VecSet(x_old[g], 0.0);
    }
    }
  // add two-grid mat and vec, if needed
  if (do_two_grid)
  {
    A[num_groups] = chi_math::PETScUtils::CreateSquareMatrix(n,N);
    chi_math::PETScUtils::InitMatrixSparsity(A[num_groups],
                                             nodal_nnz_in_diag,
                                             nodal_nnz_off_diag);
    VecDuplicate(x.front(), &x[num_groups]); // jcr needed?
//    x[num_groups] = chi_math::PETScUtils::CreateVectorWithGhosts(n,N,
//                                                        static_cast<int64_t>(ghost_dof_indices.size()),
//                                                        ghost_dof_indices);
  }

  if (do_two_grid)
    mg_diffusion::Solver::Compute_TwoGrid_VolumeFractions();

  //============================================= Create Mats and ExtVecs
  mg_diffusion::Solver::Assemble_A_bext();

  //============================================= Volume fraction for two-grid
  //                                              update

  //============================================= Field Function
  if (field_functions.empty())
  {
    for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
    {
      std::string solver_name;
      if (not TextName().empty()) solver_name = TextName() + "-";

      char buff[100];
      int dummy = snprintf(buff,4,"%03d",g);

      std::string text_name = solver_name + "Flux_g" + std::string(buff);

      using namespace chi_math;
      auto initial_field_function =
        std::make_shared<chi_physics::FieldFunction>(
          text_name,                     //Text name
          sdm_ptr,                       //Spatial Discretization
          Unknown(UnknownType::SCALAR)); //Unknown Manager

      field_functions.push_back(initial_field_function);
      chi::field_function_stack.push_back(initial_field_function);
    }//for g
  }//if not ff set

}//end initialize
