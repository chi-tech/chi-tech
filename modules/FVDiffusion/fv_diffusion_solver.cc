#include "fv_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "fv_diffusion_bndry.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

//============================================= constructor
fv_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, { {"max_iters", int64_t(500)   },
                                        {"residual_tolerance", 1.0e-2}})
{}

//============================================= destructor
fv_diffusion::Solver::~Solver()
{
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
}

//============================================= Initialize
void fv_diffusion::Solver::Initialize()
{
  chi::log.Log() << "\n"
                     << chi::program_timer.GetTimeString() << " "
                     << TextName() << ": Initializing CFEM Diffusion solver ";

  //============================================= Get grid
  grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;
  if (grid_ptr == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " No grid defined.");
 
  chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();

  uint64_t max_boundary_id = 0;
  for (const auto& id : globl_unique_bndry_ids)
    max_boundary_id = std::max(id,max_boundary_id);

  chi::log.Log() << "Max boundary id identified: " << max_boundary_id;

  for (int bndry=0; bndry<(max_boundary_id+1); bndry++)
  {
    if (boundary_preferences.find(bndry) != boundary_preferences.end())
    {
      BoundaryInfo bndry_info = boundary_preferences.at(bndry);
      auto& bndry_vals = bndry_info.second;
      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting:
        {
          boundaries.push_back({BoundaryType::Reflecting, {0.,0.,0.}});
          chi::log.Log() << "Boundary " << bndry << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1,0.0);
          boundaries.push_back({BoundaryType::Dirichlet, {bndry_vals[0],0.,0.}});
          chi::log.Log() << "Boundary " << bndry << " set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size()!=3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " Robin needs 3 values in bndry vals.");
          boundaries.push_back({BoundaryType::Robin, {bndry_vals[0],
                                                      bndry_vals[1],
                                                      bndry_vals[2]}});
          chi::log.Log() << "Boundary " << bndry << " set to robin." << bndry_vals[0]<<","<<bndry_vals[1]<<","<<bndry_vals[2];
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries.push_back({BoundaryType::Robin, {0.25,0.5,0.}});
          chi::log.Log() << "Boundary " << bndry << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size()!=3) 
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " Neumann needs 3 values in bndry vals.");
          boundaries.push_back({BoundaryType::Robin, {0.,bndry_vals[0],
                                                      bndry_vals[1]}});
          chi::log.Log() << "Boundary " << bndry << " set to neumann." << bndry_vals[0];
          break;
        }
      }//switch boundary type
    }
    else
    {
      boundaries.push_back({BoundaryType::Dirichlet, {0.,0.,0.}});
      chi::log.Log0Verbose1()
        << "No boundary preference found for boundary index " << bndry
        << "Dirichlet boundary added with zero boundary value.";
    }
  }//for bndry
  
  //============================================= Make SDM
  sdm_ptr = chi_math::SpatialDiscretization_FV::New(grid_ptr);
  const auto& sdm = *sdm_ptr;
 
  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);
 
  chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
 
  A = chi_math::PETScUtils::CreateSquareMatrix(n,N);
  x = chi_math::PETScUtils::CreateVector(n,N);
  b = chi_math::PETScUtils::CreateVector(n,N);
 
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag, OneDofPerNode);
 
  chi_math::PETScUtils::InitMatrixSparsity(A,
                                           nodal_nnz_in_diag,
                                           nodal_nnz_off_diag);

  if (field_functions_.empty())
  {
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    using namespace chi_math;
    auto initial_field_function =
      std::make_shared<chi_physics::FieldFunction>(
        text_name,                     //Text name
        sdm_ptr,                       //Spatial Discretization
        Unknown(UnknownType::SCALAR)); //Unknown/Variable

    field_functions_.push_back(initial_field_function);
    chi::field_function_stack.push_back(initial_field_function);
  }//if not ff set

}//end initialize

//========================================================== Execute
void fv_diffusion::Solver::Execute()
{
  chi::log.Log() << "\nExecuting CFEM Diffusion solver";

  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  lua_State* L = chi::console.consoleState;

  //============================================= Assemble the system
  // P ~ Present cell
  // N ~ Neighbor cell
  chi::log.Log() << "Assembling system: ";
  for (const auto& cell_P : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell_P);
    const double volume_P = cell_mapping.CellVolume(); //Volume of present cell
    const auto& x_cc_P = cell_P.centroid;

    const auto imat  = cell_P.material_id;

    const double sigma_a = CallLua_iXYZFunction(L, "Sigma_a",imat,x_cc_P);
    const double q_ext   = CallLua_iXYZFunction(L, "Q_ext"  ,imat,x_cc_P);
    const double D_P     = CallLua_iXYZFunction(L, "D_coef" ,imat,x_cc_P);

    const int64_t imap = sdm.MapDOF(cell_P, 0);
    MatSetValue(A, imap, imap, sigma_a * volume_P, ADD_VALUES);
    VecSetValue(b, imap,       q_ext   * volume_P, ADD_VALUES);

    for (size_t f=0; f < cell_P.faces.size(); ++f)
    {
      const auto& face  = cell_P.faces[f];
      const auto& x_fc  = face.centroid;
      const auto  x_PF  = x_fc - x_cc_P;
      const auto  A_f   = cell_mapping.FaceArea(f);
      const auto  A_f_n = A_f * face.normal;

      if (face.has_neighbor)
      {
        const auto& cell_N = grid.cells[face.neighbor_id];
        const int   jmat   = cell_N.material_id;
        const auto& x_cc_N = cell_N.centroid;
        const auto  x_PN   = x_cc_N - x_cc_P;

        const double D_N = CallLua_iXYZFunction(L, "D_coef",jmat,x_cc_N);

        const double w = x_PF.Norm()/x_PN.Norm();
        const double D_f = 1.0/(w/D_P + (1.0-w)/D_N);

        const double entry_ii = A_f_n.Dot(D_f * x_PN/x_PN.NormSquare());
        const double entry_ij = - entry_ii;

        const int64_t jmap = sdm.MapDOF(cell_N, 0);
        MatSetValue(A, imap, imap, entry_ii, ADD_VALUES);
        MatSetValue(A, imap, jmap, entry_ij, ADD_VALUES);
      }//internal face
      else
      {
        const auto& bndry = boundaries[face.neighbor_id];

        if (bndry.type == BoundaryType::Robin)
        {
          const auto& aval = bndry.values[0];
          const auto& bval = bndry.values[1];
          const auto& fval = bndry.values[2];

          if (std::fabs(bval) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

          if (std::fabs(aval) > 1.0e-8)
            MatSetValue(A, imap, imap, A_f * aval/bval, ADD_VALUES);
          if (std::fabs(fval) > 1.0e-8)
            VecSetValue(b, imap,       A_f * fval/bval, ADD_VALUES);
        }//if Robin

        if (bndry.type == BoundaryType::Dirichlet)
        {
          const auto& boundary_value = bndry.values[0];

          const auto& x_cc_N = x_cc_P + 2.0*x_PF;
          const auto  x_PN   = x_cc_N - x_cc_P;

          const double D_f = D_P;
          const double entry_ii = A_f_n.Dot(D_f * x_PN/x_PN.NormSquare());

          MatSetValue(A, imap, imap, entry_ii, ADD_VALUES);
          VecSetValue(b, imap,       entry_ii * boundary_value, ADD_VALUES);
        }//if Dirichlet
      }//bndry face
    }//for f
  }//for cell
 
  chi::log.Log() << "Global assembly";
 
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
 
  chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  chi::log.Log() << "Solving: ";
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
      A,               //Matrix
      TextName(),      //Solver name
      KSPCG,           //Solver type
      PCGAMG,          //Preconditioner type
      basic_options_("residual_tolerance").FloatValue(),  //Relative residual tolerance
      basic_options_("max_iters").IntegerValue()          //Max iterations
      );
 
  //============================================= Solve
  KSPSolve(petsc_solver.ksp,b,x);

  UpdateFieldFunctions();
 
  chi::log.Log() << "Done solving";

}