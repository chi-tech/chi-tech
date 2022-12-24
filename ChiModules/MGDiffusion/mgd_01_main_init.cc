#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "mg_diffusion_bndry.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//============================================= constructor
mg_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, { {"max_iters", int64_t(500)   },
                                        {"residual_tolerance", 1.0e-2}})
{}

//============================================= destructor
mg_diffusion::Solver::~Solver()
{
  for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
  {
    VecDestroy(&x[g]);
    VecDestroy(&bext[g]);
    MatDestroy(&A[g]);
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

  A.resize(mg_diffusion::Solver::num_groups, nullptr);
  bext.resize(mg_diffusion::Solver::num_groups, nullptr);
  x.resize(mg_diffusion::Solver::num_groups, nullptr);

  for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
  {
    A[g] = chi_math::PETScUtils::CreateSquareMatrix(n,N);
    x[g] = chi_math::PETScUtils::CreateVector(n,N);
    bext[g] = chi_math::PETScUtils::CreateVector(n,N);

    chi_math::PETScUtils::InitMatrixSparsity(A[g],
                                             nodal_nnz_in_diag,
                                             nodal_nnz_off_diag);
  }
  // also initalize b
  VecDuplicate(bext.front(), &b);

  //============================================= Create Mats and ExtVecs
  mg_diffusion::Solver::Assemble_A_bext();

  //============================================= Field Function
  if (field_functions.empty())
  {
    auto unk_man = OneDofPerNode;

    for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
    {
      auto initial_field_function =
        std::make_shared<chi_physics::FieldFunction>(
          std::string("mg_phi_"+std::to_string(g)),//Text name
          sdm_ptr,              //Spatial Discretization
          &x[g],                //Data vector
          unk_man);             //Unknown Manager

      field_functions.push_back(initial_field_function);
      chi::fieldfunc_stack.push_back(initial_field_function);
    }

  }//if not ff set

}//end initialize
