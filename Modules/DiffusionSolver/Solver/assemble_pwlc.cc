#include "diffusion_solver.h"

#include <PiecewiseLinear/CellViews/pwl_cellbase.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

#include <sstream>
#include "chi_log.h"

extern ChiLog chi_log;

//###################################################################
/**Assembles PWLC matrix for general cells.*/
void chi_diffusion::Solver::CFEM_Assemble_A_and_b(int cell_glob_index,
                                                  chi_mesh::Cell *cell,
                                                  int group)
{
  auto fe_view   = pwl_sdm->MapFeViewL(cell->local_id);

  //======================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id,cell,fe_view->dofs,D,q,siga,group);

  //======================================== Init cell matrix info
  typedef std::vector<double> Row;
  typedef std::vector<Row> Matrix;

  Matrix              cell_matrix;
  std::vector<double> cell_rhs;

  cell_matrix.resize(fe_view->dofs,Row(fe_view->dofs,0.0));
  cell_rhs.resize(fe_view->dofs,0.0);

  std::vector<int> dof_global_row_ind(fe_view->dofs,-1);
  std::vector<int> dof_global_col_ind(fe_view->dofs,-1);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    dof_global_row_ind[i] = pwl_sdm->MapCFEMDOF(cell->vertex_ids[i]);

    for (int j=0; j<fe_view->dofs; j++)
    {
      double mat_entry =
        D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j] +
        siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

      cell_matrix[i][j] = mat_entry;
    }//for j

    //====================== Develop RHS entry
    cell_rhs[i] = q[i]*fe_view->IntV_shapeI[i];
  }//for i
  dof_global_col_ind = dof_global_row_ind;

//  //======================================== Apply Dirichlet,Vacuum, Neumann and
//  //                                         Robin BCs
//  // Dirichlets are just collected
  std::vector<int>    dirichlet_count(fe_view->dofs,0);
  std::vector<double> dirichlet_value(fe_view->dofs,0.0);
  for (int f=0; f<cell->faces.size(); f++)
  {
    if (cell->faces[f].neighbor < 0)
    {
      int ir_boundary_index = abs(cell->faces[f].neighbor)-1;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_DIRICHLET)
      {
        auto dirichlet_bndry =
          (chi_diffusion::BoundaryDirichlet*)boundaries[ir_boundary_index];

        int num_face_dofs = cell->faces[f].vertex_ids.size();
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f][fi];
          dirichlet_count[i] += 1;
          dirichlet_value[i] += dirichlet_bndry->boundary_value;
        }
      }

      if (ir_boundary_type == DIFFUSION_ROBIN)
      {
        auto robin_bndry =
          (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];

        int num_face_dofs = cell->faces[f].vertex_ids.size();
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f][fi];

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f][fj];

            double aij = robin_bndry->a*fe_view->IntS_shapeI_shapeJ[f][i][j];
            aij /= robin_bndry->b;

            cell_matrix[i][j] += aij;
          }//for fj

          double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
          aii /= robin_bndry->b;

          cell_matrix[i][i] += aii;
        }//for fi
      }//if robin

    }//if boundary

  }//for face

  //======================================== Apply dirichlet BCs
  //Compute average dirichlet value
  for (int i=0; i<fe_view->dofs; ++i)
    dirichlet_value[i] /= (dirichlet_count[i] > 0)? dirichlet_count[i] : 1;

  for (int i=0; i<fe_view->dofs; ++i)
  {
    if (dirichlet_count[i] > 0)
    {
      cell_matrix[i].clear();
      cell_matrix[i] = std::vector<double>(fe_view->dofs,0.0);
      cell_matrix[i][i] = 1.0;
      int ir = dof_global_col_ind[i];
      MatSetValue(Aref,ir,ir,1.0,ADD_VALUES);
      dof_global_col_ind[i] = -1;
      cell_rhs[i] = dirichlet_value[i];
    }
    else
    {
      for (int j=0; j<fe_view->dofs; ++j)
      {
        if (dirichlet_count[j] > 0)
        {
          cell_rhs[i] -= cell_matrix[i][j]*dirichlet_value[j];
          cell_matrix[i][j] = 0.0;
        }
      }
    }
  }

  //======================================== Make contiguous copy of matrix
  std::vector<double> cell_matrix_cont(fe_view->dofs*fe_view->dofs,0.0);
  int n = 0;
  for (int i=0; i<fe_view->dofs; ++i)
    for (int j=0; j<fe_view->dofs; ++j)
      cell_matrix_cont[n++] = cell_matrix[i][j];

  //======================================== Add to global
  MatSetValues(Aref,
               fe_view->dofs,dof_global_row_ind.data(),
               fe_view->dofs,dof_global_col_ind.data(),
               cell_matrix_cont.data(), ADD_VALUES);

  VecSetValues(bref,
               fe_view->dofs,dof_global_row_ind.data(),
               cell_rhs.data(),ADD_VALUES);

  VecSetValues(xref,
               fe_view->dofs,
               dof_global_row_ind.data(),
               dirichlet_value.data(),INSERT_VALUES);

}