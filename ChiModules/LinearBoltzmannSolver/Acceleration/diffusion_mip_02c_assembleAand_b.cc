#include "diffusion_mip.h"
#include "acceleration.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "LinearBoltzmannSolver/lbs_structs.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"
#include "ChiConsole/chi_console.h"

#define DefaultBCDirichlet BoundaryCondition{BCType::DIRICHLET,{0,0,0}}

//###################################################################
/**Assembles both the matrix and the RHS using unit cell-matrices. These are
 * the routines used in the production versions.*/
void lbs::acceleration::DiffusionMIPSolver::
  AssembleAand_b(const std::vector<double>& q_vector)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "AssembleAand_b";
  if (m_A == nullptr or m_rhs == nullptr or m_ksp == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString() << " Starting assembly";

  const size_t num_groups   = m_uk_man.unknowns.front().num_components;

  VecSet(m_rhs, 0.0);
  for (const auto& cell : m_grid.local_cells)
  {
    const size_t num_faces    = cell.faces.size();
    const auto&  cell_mapping = m_sdm.GetCellMapping(cell);
    const size_t num_nodes    = cell_mapping.NumNodes();
    const auto   cc_nodes     = cell_mapping.GetNodeLocations();
    const auto&  unit_cell_matrices = m_unit_cell_matrices[cell.local_id];

    const auto& cell_K_matrix = unit_cell_matrices.K_matrix;
    const auto& cell_M_matrix = unit_cell_matrices.M_matrix;

    const auto& xs = m_map_mat_id_2_xs.at(cell.material_id);

    for (size_t g=0; g<num_groups; ++g)
    {
      //==================================== Get coefficient and nodal src
      const double Dg     = xs.Dg[g];
      const double sigr_g = xs.sigR[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j=0; j<num_nodes; j++)
        qg[j] = q_vector[m_sdm.MapDOFLocal(cell, j, m_uk_man, 0, g)];

      //==================================== Assemble continuous terms
      for (size_t i=0; i<num_nodes; i++)
      {
        const int64_t imap = m_sdm.MapDOF(cell,i,m_uk_man,0,g);
        double entry_rhs_i = 0.0;
        for (size_t j=0; j<num_nodes; j++)
        {
          const int64_t jmap = m_sdm.MapDOF(cell,j,m_uk_man,0,g);

          const double entry_aij =
            Dg * cell_K_matrix[i][j] +
            sigr_g * cell_M_matrix[i][j];

          entry_rhs_i += qg[j]* cell_M_matrix[i][j];

          MatSetValue(m_A, imap, jmap, entry_aij, ADD_VALUES);
        }//for j

        VecSetValue(m_rhs, imap, entry_rhs_i, ADD_VALUES);
      }//for i

      //==================================== Assemble face terms
      for (size_t f=0; f<num_faces; ++f)
      {
        const auto&  face           = cell.faces[f];
        const auto&  n_f            = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& face_M = unit_cell_matrices.face_M_matrices[f];
        const auto& face_G = unit_cell_matrices.face_G_matrices[f];
        const auto& face_Si = unit_cell_matrices.face_Si_vectors[f];

        const double hm = HPerpendicular(cell, f);

        typedef chi_mesh::MeshContinuum Grid;

        if (face.has_neighbor)
        {
          const auto&  adj_cell         = m_grid.cells[face.neighbor_id];
          const auto&  adj_cell_mapping = m_sdm.GetCellMapping(adj_cell);
          const auto   ac_nodes         = adj_cell_mapping.GetNodeLocations();
          const size_t acf              = Grid::MapCellFace(cell, adj_cell, f);
          const double hp               = HPerpendicular(adj_cell, acf);

          const auto&  adj_xs   = m_map_mat_id_2_xs.at(adj_cell.material_id);
          const double adj_Dg   = adj_xs.Dg[g];

          //========================= Compute kappa
          double kappa = 1.0;
          if (cell.Type() == chi_mesh::CellType::SLAB)
            kappa = fmax(options.penalty_factor*(adj_Dg/hp + Dg/hm)*0.5,0.25);
          if (cell.Type() == chi_mesh::CellType::POLYGON)
            kappa = fmax(options.penalty_factor*(adj_Dg/hp + Dg/hm)*0.5,0.25);
          if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
            kappa = fmax(options.penalty_factor*2.0*(adj_Dg/hp + Dg/hm)*0.5,0.25);

          //========================= Assembly penalty terms
          for (size_t fi=0; fi<num_face_nodes; ++fi)
          {
            const int i  = cell_mapping.MapFaceNode(f,fi);
            const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

            for (size_t fj=0; fj<num_face_nodes; ++fj)
            {
              const int jm = cell_mapping.MapFaceNode(f,fj);      //j-minus
              const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                             f, acf, fj);         //j-plus
              const int64_t jmmap = m_sdm.MapDOF(cell, jm, m_uk_man, 0, g);
              const int64_t jpmap = m_sdm.MapDOF(adj_cell, jp, m_uk_man, 0, g);

              const double aij = kappa * face_M[i][jm];

              MatSetValue(m_A, imap, jmmap, aij,ADD_VALUES);
              MatSetValue(m_A, imap, jpmap,-aij,ADD_VALUES);
            }//for fj
          }//for fi

          //========================= Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (int i=0; i<num_nodes; i++)
          {
            const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

            for (int fj=0; fj<num_face_nodes; fj++)
            {
              const int jm = cell_mapping.MapFaceNode(f,fj);      //j-minus
              const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                             f, acf, fj);         //j-plus
              const int64_t jmmap = m_sdm.MapDOF(cell,jm,m_uk_man,0,g);
              const int64_t jpmap = m_sdm.MapDOF(adj_cell,jp,m_uk_man,0,g);

              const double aij = -0.5*Dg*n_f.Dot(face_G[jm][i]);

              MatSetValue(m_A, imap, jmmap, aij, ADD_VALUES);
              MatSetValue(m_A, imap, jpmap,-aij, ADD_VALUES);
            }//for fj
          }//for i

          // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
          for (int fi=0; fi<num_face_nodes; fi++)
          {
            const int im = cell_mapping.MapFaceNode(f,fi);       //i-minus
            const int ip = MapFaceNodeDisc(cell,adj_cell,cc_nodes,ac_nodes,
                                           f,acf,fi);            //i-plus
            const int64_t immap = m_sdm.MapDOF(cell,im,m_uk_man,0,g);
            const int64_t ipmap = m_sdm.MapDOF(adj_cell,ip,m_uk_man,0,g);

            for (int j=0; j<num_nodes; j++)
            {
              const int64_t jmap = m_sdm.MapDOF(cell, j, m_uk_man, 0, g);

              const double aij = -0.5*Dg*n_f.Dot(face_G[im][j]);

              MatSetValue(m_A,immap,jmap, aij,ADD_VALUES);
              MatSetValue(m_A,ipmap,jmap,-aij,ADD_VALUES);
            }//for j
          }//for fi

        }//internal face
        else
        {
          auto bc = DefaultBCDirichlet;
          try {bc = m_bcs.at(face.neighbor_id);}
          catch (const std::out_of_range& oor)
          {throw std::logic_error(fname + ": unmapped boundary id.");}

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            //========================= Compute kappa
            double kappa = 1.0;
            if (cell.Type() == chi_mesh::CellType::SLAB)
              kappa = fmax(options.penalty_factor*Dg/hm,0.25);
            if (cell.Type() == chi_mesh::CellType::POLYGON)
              kappa = fmax(options.penalty_factor*Dg/hm,0.25);
            if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor*2.0*Dg/hm,0.25);

            //========================= Assembly penalty terms
            for (size_t fi=0; fi<num_face_nodes; ++fi)
            {
              const int i  = cell_mapping.MapFaceNode(f,fi);
              const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              for (size_t fj=0; fj<num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f,fj);
                const int64_t jmmap = m_sdm.MapDOF(cell, jm, m_uk_man, 0, g);

                const double aij = kappa*face_M[i][jm];
                const double aij_bc_value = aij*bc_value;

                MatSetValue(m_A  , imap, jmmap, aij,ADD_VALUES);
                VecSetValue(m_rhs, imap, aij_bc_value, ADD_VALUES);
              }//for fj
            }//for fi

            //========================= Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = n dot nabla bk

            // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i=0; i<num_nodes; i++)
            {
              const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              for (size_t j=0; j<num_nodes; j++)
              {
                const int64_t jmap = m_sdm.MapDOF(cell, j, m_uk_man, 0, g);

                const double aij = -Dg*n_f.Dot(face_G[j][i] + face_G[i][j]);
                const double aij_bc_value = aij*bc_value;


                MatSetValue(m_A, imap, jmap, aij, ADD_VALUES);
                VecSetValue(m_rhs, imap, aij_bc_value, ADD_VALUES);
              }//for fj
            }//for i
          }//Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double aval = bc.values[0];
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12) continue; //a and f assumed zero

            for (size_t fi=0; fi<num_face_nodes; fi++)
            {
              const int i  = cell_mapping.MapFaceNode(f,fi);
              const int64_t ir = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              if (std::fabs(aval) >= 1.0e-12)
              {
                for (size_t fj=0; fj<num_face_nodes; fj++)
                {
                  const int j  = cell_mapping.MapFaceNode(f,fj);
                  const int64_t jr = m_sdm.MapDOF(cell, j, m_uk_man, 0, g);

                  const double aij = (aval/bval) * face_M[i][j];

                  MatSetValue(m_A,ir ,jr, aij,ADD_VALUES);
                }//for fj
              }//if a nonzero

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval/bval) * face_Si[i];

                VecSetValue(m_rhs,ir, rhs_val, ADD_VALUES);
              }//if f nonzero
            }//for fi
          }//Robin BC
        }//boundary face
      }//for face
    }//for g
  }//for cell

  MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(m_rhs);
  VecAssemblyEnd(m_rhs);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(m_A, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  KSPSetOperators(m_ksp, m_A, m_A);

  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString() << " Assembly completed";

  PC pc;
  KSPGetPC(m_ksp, &pc);
  PCSetUp(pc);

  KSPSetUp(m_ksp);
}