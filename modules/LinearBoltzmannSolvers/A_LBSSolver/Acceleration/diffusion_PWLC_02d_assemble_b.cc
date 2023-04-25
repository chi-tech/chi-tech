#include "diffusion_PWLC.h"
#include "acceleration.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "ChiPhysics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"
#include "ChiConsole/chi_console.h"

#define DefaultBCDirichlet                                                     \
  BoundaryCondition                                                            \
  {                                                                            \
    BCType::DIRICHLET, { 0, 0, 0 }                                             \
  }

// ###################################################################
/**Assembles the RHS using unit cell-matrices. These are
 * the routines used in the production versions.*/
void lbs::acceleration::DiffusionPWLCSolver::Assemble_b(
  const std::vector<double>& q_vector)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "Assemble_b";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString()
                   << " Starting assembly";

  const size_t num_groups = uk_man_.unknowns_.front().num_components_;

  VecSet(rhs_, 0.0);
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto& unit_cell_matrices = unit_cell_matrices_[cell.local_id_];

    const auto& cell_M_matrix = unit_cell_matrices.M_matrix;
    const auto& cell_Vi = unit_cell_matrices.Vi_vectors;

    //=========================================== Mark dirichlet nodes
    std::vector<bool> node_is_dirichlet(num_nodes, false);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (not face.has_neighbor_)
      {
        auto bc = DefaultBCDirichlet;
        if (bcs_.count(face.neighbor_id_) > 0) bc = bcs_.at(face.neighbor_id_);

        if (bc.type != BCType::DIRICHLET) continue;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
          node_is_dirichlet[cell_mapping.MapFaceNode(f, fi)] = true;
      }
    }

    for (size_t g = 0; g < num_groups; ++g)
    {
      //==================================== Get coefficient and nodal src
      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      //==================================== Assemble continuous terms
      for (size_t i = 0; i < num_nodes; i++)
      {
        if (node_is_dirichlet[i]) continue;
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; j++)
          entry_rhs_i += qg[j] * cell_M_matrix[i][j];

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      //==================================== Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces_[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& face_Si = unit_cell_matrices.face_Si_vectors[f];

        if (not face.has_neighbor_)
        {
          auto bc = DefaultBCDirichlet;
          if (bcs_.count(face.neighbor_id_) > 0)
            bc = bcs_.at(face.neighbor_id_);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            //========================= Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              VecSetValue(rhs_, imap, bc_value * cell_Vi[i], ADD_VALUES);
            } // for fi

          } // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12) continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval / bval) * face_Si[i];

                VecSetValue(rhs_, ir, rhs_val, ADD_VALUES);
              } // if f nonzero
            }   // for fi
          }     // Robin BC
        }       // boundary face
      }         // for face
    }           // for g
  }             // for cell

  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);

  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString()
                   << " Assembly completed";
}

// ###################################################################
/**Assembles the RHS using unit cell-matrices. These are
 * the routines used in the production versions.*/
void lbs::acceleration::DiffusionPWLCSolver::Assemble_b(Vec petsc_q_vector)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "Assemble_b";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString()
                   << " Starting assembly";

  const size_t num_groups = uk_man_.unknowns_.front().num_components_;

  const double* q_vector;
  VecGetArrayRead(petsc_q_vector, &q_vector);

  VecSet(rhs_, 0.0);
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto& unit_cell_matrices = unit_cell_matrices_[cell.local_id_];

    const auto& cell_M_matrix = unit_cell_matrices.M_matrix;
    const auto& cell_Vi = unit_cell_matrices.Vi_vectors;

    //=========================================== Mark dirichlet nodes
    std::vector<bool> node_is_dirichlet(num_nodes, false);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (not face.has_neighbor_)
      {
        auto bc = DefaultBCDirichlet;
        if (bcs_.count(face.neighbor_id_) > 0) bc = bcs_.at(face.neighbor_id_);

        if (bc.type != BCType::DIRICHLET) continue;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
          node_is_dirichlet[cell_mapping.MapFaceNode(f, fi)] = true;
      }
    }

    for (size_t g = 0; g < num_groups; ++g)
    {
      //==================================== Get coefficient and nodal src
      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      //==================================== Assemble continuous terms
      for (size_t i = 0; i < num_nodes; i++)
      {
        if (node_is_dirichlet[i]) continue;
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; j++)
          entry_rhs_i += qg[j] * cell_M_matrix[i][j];

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      //==================================== Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces_[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& face_Si = unit_cell_matrices.face_Si_vectors[f];

        if (not face.has_neighbor_)
        {
          auto bc = DefaultBCDirichlet;
          if (bcs_.count(face.neighbor_id_) > 0)
            bc = bcs_.at(face.neighbor_id_);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            //========================= Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              VecSetValue(rhs_, imap, bc_value * cell_Vi[i], ADD_VALUES);
            } // for fi

          } // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12) continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval / bval) * face_Si[i];

                VecSetValue(rhs_, ir, rhs_val, ADD_VALUES);
              } // if f nonzero
            }   // for fi
          }     // Robin BC
        }       // boundary face
      }         // for face
    }           // for g
  }             // for cell

  VecRestoreArrayRead(petsc_q_vector, &q_vector);

  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);

  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString()
                   << " Assembly completed";
}