#include "diffusion_mip.h"

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

void lbs::acceleration::DiffusionMIPSolver::
  Assemble_b2(const std::vector<double>& q_vector)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "Assemble_b";
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

    const auto& cell_M_matrix = unit_cell_matrices.M_matrix;

    const auto& xs = m_map_mat_id_2_xs.at(cell.material_id);

    for (size_t g=0; g<num_groups; ++g)
    {
      //==================================== Get coefficient and nodal src
      const double Dg     = xs.Dg[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j=0; j<num_nodes; j++)
        qg[j] = q_vector[m_sdm.MapDOFLocal(cell, j, m_uk_man, 0, g)];

      //==================================== Assemble continuous terms
      for (size_t i=0; i<num_nodes; i++)
      {
        const int64_t imap = m_sdm.MapDOF(cell,i,m_uk_man,0,g);
        double entry_rhs_i = 0.0;
        for (size_t j=0; j<num_nodes; j++)
          entry_rhs_i += qg[j]* cell_M_matrix[i][j];

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

        if (not face.has_neighbor)
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
              kappa = fmax(4.0*Dg/hm,0.25);
            if (cell.Type() == chi_mesh::CellType::POLYGON)
              kappa = fmax(4.0*Dg/hm,0.25);
            if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
              kappa = fmax(8.0*Dg/hm,0.25);

            //========================= Assembly penalty terms
            for (size_t fi=0; fi<num_face_nodes; ++fi)
            {
              const int i  = cell_mapping.MapFaceNode(f,fi);
              const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              for (size_t fj=0; fj<num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f,fj);

                const double aij = kappa*face_M[i][jm];
                const double aij_bc_value = aij*bc_value;

                VecSetValue(m_rhs, imap, aij_bc_value, ADD_VALUES);
              }//for fj
            }//for fi

            //========================= Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = 0.5* n dot nabla bk

            // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i=0; i<num_nodes; i++)
            {
              const int64_t imap = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              for (size_t j=0; j<num_nodes; j++)
              {
                const double aij = -0.5*Dg*n_f.Dot(face_G[j][i] + face_G[i][j]);
                const double aij_bc_value = aij*bc_value;

                VecSetValue(m_rhs, imap, aij_bc_value, ADD_VALUES);
              }//for fj
            }//for i
          }//Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12) continue; //a and f assumed zero

            for (size_t fi=0; fi<num_face_nodes; fi++)
            {
              const int i  = cell_mapping.MapFaceNode(f,fi);
              const int64_t ir = m_sdm.MapDOF(cell, i, m_uk_man, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval/bval) * face_Si[i];

                VecSetValue(m_rhs,ir, -rhs_val, ADD_VALUES);
              }//if f nonzero
            }//for fi
          }//Robin BC
        }//boundary face
      }//for face
    }//for g
  }//for cell

  VecAssemblyBegin(m_rhs);
  VecAssemblyEnd(m_rhs);

  if (options.verbose)
    chi::log.Log() << chi::program_timer.GetTimeString() << " Assembly completed";
}