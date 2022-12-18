#include "dfem_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "dfem_diffusion_bndry.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

// jcr: which ones? mesh? meshhandler? meshcontinuum?unk_man?
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/chi_mesh.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#define DefaultBCDirichlet BoundaryCondition{BCType::DIRICHLET,{0,0,0}}

//============================================= constructor
dfem_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, { {"max_iters", int64_t(500)   },
                                        {"residual_tolerance", 1.0e-2}}
                   )
{}

//============================================= destructor
dfem_diffusion::Solver::~Solver()
{
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
}

//============================================= Initialize
void dfem_diffusion::Solver::Initialize(bool verbose)
{
  chi::log.Log() << "\n"
                 << chi::program_timer.GetTimeString() << " "
                 << TextName() << ": Initializing DFEM Diffusion solver ";
  this->verbose_info = verbose;
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
          chi::log.Log() << "Boundary " << bndry << " set to robin."
                         << bndry_vals[0]<<","<<bndry_vals[1]<<","<<bndry_vals[2];
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
  sdm_ptr = chi_math::SpatialDiscretization_PWLD::New(grid_ptr);
  const auto& sdm = *sdm_ptr;
 
  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  // jcr: do I need static_cast<int64_t>?
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
  if (field_functions.empty())
  {
    auto unk_man = OneDofPerNode;
    field.resize(n);
    auto initial_field_function =
      std::make_shared<chi_physics::FieldFunction>(
        std::string("phi"),   //Text name
        sdm_ptr,              //Spatial Discretization
        &field,               //Data vector
        unk_man);             //Unknown Manager

      field_functions.push_back(initial_field_function);
      chi::fieldfunc_stack.push_back(initial_field_function);
  }//if not ff set

}//end initialize

//========================================================== Execute
void dfem_diffusion::Solver::Execute()
{
  chi::log.Log() << "\nExecuting DFEM IP Diffusion solver";

  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  lua_State* L = chi::console.consoleState;

  //============================================= Assemble the system
  // is this needed?
  VecSet(b, 0.0);

  chi::log.Log() << "Assembling system: ";

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes   = cell_mapping.NumNodes();
    const auto   cc_nodes    = cell_mapping.GetNodeLocations();
    const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();

//    cout << std::endl;
//    chi::log.Log() << "cell/node layout : " << icell;
//    for (int i = 0; i < cc_nodes.size(); ++i){
//      for (int j = 0; j < 3; ++j) { // jcr: hardcoded 3, ask Jan about .size for vec3
//        cout << cc_nodes[i][j] << ", ";
//      }
//      cout << std::endl;
//    }

    const auto imat  = cell.material_id;
    MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl cell_rhs(num_nodes, 0.0);

    //==================================== Assemble volumetric terms
    for (size_t i=0; i<num_nodes; ++i)
    {
      size_t imap = sdm.MapDOF(cell, i);

      for (size_t j=0; j<num_nodes; ++j)
      {
        size_t jmap = sdm.MapDOF(cell, j);
        double entry_aij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_aij +=
            (
              CallLua_iXYZFunction(L, "D_coef",imat,qp_data.QPointXYZ(qp))  *
              qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp))
              +
              CallLua_iXYZFunction(L, "Sigma_a",imat,qp_data.QPointXYZ(qp))  *
              qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp)
            )
            *
            qp_data.JxW(qp);
        }//for qp
        MatSetValue(A, imap, jmap, entry_aij, ADD_VALUES);
      }//for j
      double entry_rhs_i = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        entry_rhs_i += CallLua_iXYZFunction(L,"Q_ext",imat,qp_data.QPointXYZ(qp))
          * qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      VecSetValue(b, imap, entry_rhs_i, ADD_VALUES);
    }//for i

    //==================================== Assemble face terms
    const size_t num_faces = cell.faces.size();
    for (size_t f = 0; f < num_faces; ++f) {
      const auto &face = cell.faces[f];
      const auto &n_f = face.normal;
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      const auto fqp_data = cell_mapping.MakeFaceQuadraturePointData(f);

      const double hm = HPerpendicular(cell, f);

      typedef chi_mesh::MeshContinuum Grid;

      // interior face
      if (face.has_neighbor) {
        const auto &adj_cell = grid.cells[face.neighbor_id];
        const auto &adj_cell_mapping = sdm.GetCellMapping(adj_cell);
        const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
        const size_t acf = Grid::MapCellFace(cell, adj_cell, f);
        const double hp_neigh = HPerpendicular(adj_cell, acf);

        const auto imat_neigh = adj_cell.material_id;

        //========================= Compute Ckappa IP
        double Ckappa = 1.0;
        if (cell.Type() == chi_mesh::CellType::SLAB)
          Ckappa = 2.0;
        if (cell.Type() == chi_mesh::CellType::POLYGON)
          Ckappa = 2.0;
        if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
          Ckappa = 4.0;

        //========================= Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi) {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const size_t imap = sdm.MapDOF(cell, i);

          for (size_t fj = 0; fj < num_face_nodes; ++fj) {
            const int jm = cell_mapping.MapFaceNode(f, fj);      //j-minus
            const size_t jmmap = sdm.MapDOF(cell, jm);

            const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                           f, acf, fj);         //j-plus
            const size_t jpmap = sdm.MapDOF(adj_cell, jp);

            double aij = 0.0;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              aij += Ckappa *
                     (CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) / hm +
                      CallLua_iXYZFunction(L, "D_coef", imat_neigh, fqp_data.QPointXYZ(qp)) / hp_neigh) / 2.
                     *
                     fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) *
                     fqp_data.JxW(qp);

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          }//for fj
        }//for fi

        //========================= Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // {{D d_n b_i}}[[Phi]]
        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-

        // loop over node of current cell (gradient of b_i)
        for (int i = 0; i < num_nodes; ++i) {
          const size_t imap = sdm.MapDOF(cell, i);

          // loop over faces
          for (int fj = 0; fj < num_face_nodes; ++fj) {
            const int jm = cell_mapping.MapFaceNode(f, fj);      //j-minus
            const size_t jmmap = sdm.MapDOF(cell, jm);
            const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                           f, acf, fj);         //j-plus
            const size_t jpmap = sdm.MapDOF(adj_cell, jp);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              vec_aij +=
                CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) *
                fqp_data.ShapeValue(jm, qp) * fqp_data.ShapeGrad(i, qp) *
                fqp_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij);  // jcr where is this minus from??

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          }//for fj
        }//for i

        // {{D d_n Phi}}[[b_i]]
        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
        for (int fi = 0; fi < num_face_nodes; ++fi) {
          const int im = cell_mapping.MapFaceNode(f, fi);       //i-minus
          const size_t immap = sdm.MapDOF(cell, im);

          const int ip = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                         f, acf, fi);            //i-plus
          const size_t ipmap = sdm.MapDOF(adj_cell, ip);

          for (int j = 0; j < num_nodes; ++j) {
            const size_t jmap = sdm.MapDOF(cell, j);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              vec_aij +=
                CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) *
                fqp_data.ShapeValue(im, qp) * fqp_data.ShapeGrad(j, qp) *
                fqp_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij); // jcr where is this minus from??

            MatSetValue(A, immap, jmap, aij, ADD_VALUES);
            MatSetValue(A, ipmap, jmap, -aij, ADD_VALUES);
          }//for j
        }//for fi

      }//internal face
      else { // boundary face
        const auto &bndry = boundaries[face.neighbor_id];
        // Robin boundary
        if (bndry.type == BoundaryType::Robin) {
          const auto qp_face_data = cell_mapping.MakeFaceQuadraturePointData(f);
          const size_t num_face_nodes = face.vertex_ids.size();

          const auto &aval = bndry.values[0];
          const auto &bval = bndry.values[1];
          const auto &fval = bndry.values[2];

          chi::log.Log() << "Boundary  set as Robin with a,b,f = ("
                         << aval << ","
                         << bval << ","
                         << fval << ") ";
          // true Robin when a!=0, otherwise, it is a Neumann:
          // Assert if b=0
          if (std::fabs(bval) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

          for (size_t fi = 0; fi < num_face_nodes; fi++) {
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const size_t ir = sdm.MapDOF(cell, i);

            if (std::fabs(aval) >= 1.0e-12) {
              for (size_t fj = 0; fj < num_face_nodes; fj++) {
                const uint j = cell_mapping.MapFaceNode(f, fj);
                const size_t jr = sdm.MapDOF(cell, j);

                double aij = 0.0;
                for (size_t qp: fqp_data.QuadraturePointIndices())
                  aij += fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(j, qp) *
                         fqp_data.JxW(qp);
                aij *= (aval / bval);

                MatSetValue(A, ir, jr, aij, ADD_VALUES);
              }//for fj
            }//if a nonzero

            if (std::fabs(fval) >= 1.0e-12) {
              double rhs_val = 0.0;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                rhs_val += fqp_data.ShapeValue(i, qp) * fqp_data.JxW(qp);
              rhs_val *= (fval / bval);

              VecSetValue(b, ir, rhs_val, ADD_VALUES);
            }//if f nonzero
          }//for fi
        }//Robin BC
        else if (bndry.type == BoundaryType::Dirichlet) {
          const double bc_value = bndry.values[0];
          //========================= Compute kappa
          double Ckappa = 2.0;
          if (cell.Type() == chi_mesh::CellType::SLAB)
            Ckappa = 4.0; // fmax(4.0*Dg/hm,0.25);
          if (cell.Type() == chi_mesh::CellType::POLYGON)
            Ckappa = 4.0;
          if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
            Ckappa = 8.0;

          //========================= Assembly penalty terms
          for (size_t fi = 0; fi < num_face_nodes; ++fi) {
            // jcr int, uint, size_t ?
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const size_t imap = sdm.MapDOF(cell, i);

            for (size_t fj = 0; fj < num_face_nodes; ++fj) {
              const uint jm = cell_mapping.MapFaceNode(f, fj);
              const size_t jmmap = sdm.MapDOF(cell, jm);

              double aij = 0.0;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                aij += Ckappa *
                       CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) / hm *
                       fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) *
                       fqp_data.JxW(qp);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
              VecSetValue(b, imap, aij_bc_value, ADD_VALUES);
            }//for fj
          }//for fi

          //========================= Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (size_t i = 0; i < num_nodes; i++) {
            // jcr: ask size_t versus int64_t
            const size_t imap = sdm.MapDOF(cell, i);

            for (size_t j = 0; j < num_nodes; j++) {
              const size_t jmap = sdm.MapDOF(cell, j);

              chi_mesh::Vector3 vec_aij;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                vec_aij +=
                  (fqp_data.ShapeValue(j, qp) * fqp_data.ShapeGrad(i, qp) +
                   fqp_data.ShapeValue(i, qp) * fqp_data.ShapeGrad(j, qp)) *
                  fqp_data.JxW(qp) *
                  CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp));

              const double aij = -n_f.Dot(vec_aij);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A, imap, jmap, aij, ADD_VALUES);
              VecSetValue(b, imap, aij_bc_value, ADD_VALUES);
            }//for fj
          }//for i
        }//Dirichlet BC
        else {

        }//else BC
      }//boundary face
    }//for face f
  }//for cell

 
  chi::log.Log() << "Global assembly";
 
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

//  MatView(A, PETSC_VIEWER_STDERR_WORLD);
//
//  PetscViewer viewer;
//  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A.m",&viewer);
//  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//  MatView(A,viewer);
//  PetscViewerPopFormat(viewer);
//  PetscViewerDestroy(&viewer);

  chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  chi::log.Log() << "Solving: ";
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
      A,               //Matrix
      TextName(),      //Solver name
      KSPCG,           //Solver type
      PCGAMG,          //Preconditioner type
      basic_options("residual_tolerance").FloatValue(),  //Relative residual tolerance
      basic_options("max_iters").IntegerValue()          //Max iterations
      );
 
  //============================================= Solve
  KSPSolve(petsc_solver.ksp,b,x);
 
  chi::log.Log() << "Done solving";

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  sdm.LocalizePETScVector(x,field,OneDofPerNode);


}