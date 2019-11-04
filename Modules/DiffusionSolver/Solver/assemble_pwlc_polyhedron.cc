#include "diffusion_solver.h"

#include <PiecewiseLinear/CellViews/pwl_polyhedron.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::CFEM_Ab_Polyhedron(int cell_glob_index,
                                               chi_mesh::Cell *cell,
                                               int group)
{
  chi_mesh::CellPolyhedron* polyh_cell =
    (chi_mesh::CellPolyhedron*)(cell);
  PolyhedronFEView* fe_view =
    (PolyhedronFEView*)pwl_discr->MapFeView(cell_glob_index);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,group);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = mesher->MapNode(polyh_cell->v_indices[i]);

    int ir_boundary_type;
    if (!ApplyDirichletI(ir,&ir_boundary_type))
    {
      //====================== Develop matrix entry
      for (int j=0; j<fe_view->dofs; j++)
      {
        int jr =  mesher->MapNode(polyh_cell->v_indices[j]);
        double jr_mat_entry =
          D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

        jr_mat_entry +=
          siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

        int jr_boundary_type;
        if (!ApplyDirichletJ(jr,ir,jr_mat_entry,&jr_boundary_type))
        {
          MatSetValue(Aref,ir,jr,jr_mat_entry,ADD_VALUES);
        }
      }//for j

      //====================== Develop RHS entry
      double rhsvalue =0.0;
      rhsvalue = q[i]*fe_view->IntV_shapeI[i];
      VecSetValue(bref,ir,rhsvalue,ADD_VALUES);
    }//if ir not dirichlet

  }//for i

  //======================================== Apply Vacuum, Neumann and Robin
  //                                         BCs
  for (int f=0; f<polyh_cell->faces.size(); f++)
  {
    int num_face_dofs = polyh_cell->faces[f]->v_indices.size();
    if (polyh_cell->faces[f]->face_indices[NEIGHBOR] < 0)
    {
      int ir_boundary_index =
        abs(polyh_cell->faces[f]->face_indices[NEIGHBOR])-1;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_ROBIN)
      {
        auto robin_bndry =
          (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];

        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f]->cell_dof[fi];
          int ir = mesher->MapNode(polyh_cell->v_indices[i]);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f]->cell_dof[fj];
            int jr = mesher->MapNode(polyh_cell->v_indices[j]);

            double aij = robin_bndry->a*fe_view->IntS_shapeI_shapeJ[f][i][j];
            aij /= robin_bndry->b;

            MatSetValue(Aref,ir ,jr, aij,ADD_VALUES);
          }//for fj

          double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
          aii /= robin_bndry->b;

          MatSetValue(Aref,ir ,ir, aii,ADD_VALUES);
        }//for fi
      }//if vacuum
    }//if boundary

  }//for face
}

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::CFEM_Assemble_A_and_b(int cell_glob_index,
                                                  chi_mesh::Cell *cell,
                                                  int group)
{
  chi_mesh::CellBase* base_cell =
    (chi_mesh::CellBase*)(cell);
  PolyhedronFEView* fe_view =
    (PolyhedronFEView*)pwl_discr->MapFeView(cell_glob_index);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,group);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = mesher->MapNode(base_cell->vertex_ids[i]);

    int ir_boundary_type;
    if (!ApplyDirichletI(ir,&ir_boundary_type))
    {
      //====================== Develop matrix entry
      for (int j=0; j<fe_view->dofs; j++)
      {
        int jr =  mesher->MapNode(base_cell->vertex_ids[j]);
        double jr_mat_entry =
          D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

        jr_mat_entry +=
          siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

        int jr_boundary_type;
        if (!ApplyDirichletJ(jr,ir,jr_mat_entry,&jr_boundary_type))
        {
          MatSetValue(Aref,ir,jr,jr_mat_entry,ADD_VALUES);
        }
      }//for j

      //====================== Develop RHS entry
      double rhsvalue =0.0;
      rhsvalue = q[i]*fe_view->IntV_shapeI[i];
      VecSetValue(bref,ir,rhsvalue,ADD_VALUES);
    }//if ir not dirichlet

  }//for i

  //======================================== Apply Vacuum, Neumann and Robin
  //                                         BCs
  for (int f=0; f<base_cell->faces.size(); f++)
  {
    int num_face_dofs = base_cell->faces[f].vertex_ids.size();
    if (base_cell->faces[f].neighbor < 0)
    {
      int ir_boundary_index =
        abs(base_cell->faces[f].neighbor)-1;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_ROBIN)
      {
        auto robin_bndry =
          (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];

        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f]->cell_dof[fi];
          int ir = mesher->MapNode(base_cell->vertex_ids[i]);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f]->cell_dof[fj];
            int jr = mesher->MapNode(base_cell->vertex_ids[j]);

            double aij = robin_bndry->a*fe_view->IntS_shapeI_shapeJ[f][i][j];
            aij /= robin_bndry->b;

            MatSetValue(Aref,ir ,jr, aij,ADD_VALUES);
          }//for fj

          double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
          aii /= robin_bndry->b;

          MatSetValue(Aref,ir ,ir, aii,ADD_VALUES);
        }//for fi
      }//if vacuum
    }//if boundary

  }//for face
}