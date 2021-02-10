#include "diffusion_solver.h"

#include "ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_cellbase.h"

#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_reflecting.h"
#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_robin.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::PWLD_Assemble_A_and_b_GAGG(
                                               int cell_glob_index,
                                               chi_mesh::Cell *cell,
                                               DiffusionIPCellView* cell_ip_view)
{
  auto fe_view = (CellPWLFEView*)pwl_sdm->MapFeViewL(cell->local_id);

  for (int gr=0; gr<G; gr++)
  {
    //====================================== Process material id
    int mat_id = cell->material_id;

    std::vector<double> D(fe_view->dofs,1.0);
    std::vector<double> q(fe_view->dofs,1.0);
    std::vector<double> siga(fe_view->dofs,1.0);

    GetMaterialProperties(mat_id,cell,fe_view->dofs,D,q,siga,gi+gr);

    //========================================= Loop over DOFs
    for (int i=0; i<fe_view->dofs; i++)
    {
      int ir = cell_ip_view->MapDof(i);
      int ig = cell->vertex_ids[i];
      double rhsvalue =0.0;

      int ir_boundary_type;
      //if (!ApplyDirichletI(ir*G+gr,&ir_boundary_type,ig))
      {
        //====================== Develop matrix entry
        for (int j=0; j<fe_view->dofs; j++)
        {
          int jr =  cell_ip_view->MapDof(j);
          int jg = cell->vertex_ids[j];
          double jr_mat_entry =
            D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

          jr_mat_entry +=
            siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

          int jr_boundary_type;
          //if (!ApplyDirichletJ(jr*G+gr,ir*G+gr,jr_mat_entry,&jr_boundary_type,jg))
          {
            MatSetValue(Aref,ir*G+gr,jr*G+gr,jr_mat_entry,ADD_VALUES);
          }

          rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
        }//for j

        //====================== Apply RHS entry
        VecSetValue(bref,ir*G+gr,rhsvalue,ADD_VALUES);
      }//if ir not dirichlet

    }//for i


    //========================================= Loop over faces
    int num_faces = cell->faces.size();
    for (int f=0; f<num_faces; f++)
    {
      auto& face = cell->faces[f];
      int neighbor = face.neighbor;

      //================================== Get face normal
      chi_mesh::Vector3 n  = face.normal;

      int num_face_dofs = face.vertex_ids.size();

      if (neighbor >=0)
      {
        chi_mesh::Cell*           adj_cell    = nullptr;
        CellPWLFEView*               adj_fe_view = nullptr;
        DiffusionIPCellView*     adj_ip_view  = nullptr;
        int                              fmap = -1;

        //========================= Get adj cell information
        if (cell->faces[f].IsNeighborLocal(grid))  //Local
        {
          int adj_cell_local_index = face.GetNeighborLocalID(grid);
          adj_cell      = &grid->local_cells[adj_cell_local_index];
          adj_ip_view   = ip_cell_views[adj_cell_local_index];
          adj_fe_view   = (CellPWLFEView*)pwl_sdm->MapFeViewL(adj_cell_local_index);
        }//local
        else //Non-local
        {
          int locI = face.GetNeighborPartitionID(grid);
          adj_ip_view = GetBorderIPView(locI,neighbor);
          adj_cell    = (chi_mesh::Cell*)GetBorderCell(locI,neighbor);
          adj_fe_view = (CellPWLFEView*)GetBorderFEView(locI,neighbor);
        }//non-local

        //========================= Check valid information
        if (adj_cell == nullptr || adj_fe_view == nullptr ||
            adj_ip_view == nullptr)
        {
          chi_log.Log(LOG_ALL)
            << "Error in MIP cell information.";
          exit(EXIT_FAILURE);
        }

        //========================= Get the current map to the adj cell's face
        fmap = MapCellFace(cell, adj_cell, f);

        //========================= Compute penalty coefficient
        double hp = HPerpendicular(adj_cell, adj_fe_view, fmap);
        double hm = HPerpendicular(cell, fe_view, f);

        std::vector<double> adj_D,adj_Q,adj_sigma;

        GetMaterialProperties(adj_cell->material_id,
                              adj_cell,
                              adj_fe_view->dofs,
                              adj_D,
                              adj_Q,
                              adj_sigma,
                              gi+gr);

        //========================= Compute surface average D
        double D_avg = 0.0;
        double intS = 0.0;
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i = fe_view->face_dof_mappings[f][fi];
          D_avg += D[i]*fe_view->IntS_shapeI[i][f];
          intS += fe_view->IntS_shapeI[i][f];
        }
        D_avg /= intS;

        //========================= Compute surface average D_adj
        double adj_D_avg = 0.0;
        double adj_intS = 0.0;
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i    = fe_view->face_dof_mappings[f][fi];
          int imap = MapCellDof(adj_cell, cell->vertex_ids[i]);
          adj_D_avg += adj_D[imap]*adj_fe_view->IntS_shapeI[imap][fmap];
          adj_intS += adj_fe_view->IntS_shapeI[imap][fmap];
        }
        adj_D_avg /= adj_intS;

        //========================= Compute kappa
        double kappa = 1.0;
        if (cell->Type() == chi_mesh::CellType::SLAB)
          kappa = fmax(2.0*(adj_D_avg/hp + D_avg/hm),0.25);
        if (cell->Type() == chi_mesh::CellType::POLYGON)
          kappa = fmax(2.0*(adj_D_avg/hp + D_avg/hm),0.25);
        if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
          kappa = fmax(4.0*(adj_D_avg/hp + D_avg/hm),0.25);

        //========================= Assembly penalty terms
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f][fi];
          int ir = cell_ip_view->MapDof(i);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f][fj];
            int jr = cell_ip_view->MapDof(j);
            int jmap  = MapCellDof(adj_cell, cell->faces[f].vertex_ids[fj]);
            int jrmap = adj_ip_view->MapDof(jmap);

            double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

            MatSetValue(Aref,ir*G+gr,jr   *G+gr, aij,ADD_VALUES);
            MatSetValue(Aref,ir*G+gr,jrmap*G+gr,-aij,ADD_VALUES);
          }//for fj

        }//for fi


        //========================= Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // -Di^- bj^- and
        // -Dj^- bi^-
        for (int i=0; i<fe_view->dofs; i++)
        {
          int ir = cell_ip_view->MapDof(i);

          for (int j=0; j<fe_view->dofs; j++)
          {
            int jr = cell_ip_view->MapDof(j);

            double gij =
              n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][i][j] +
                    fe_view->IntS_shapeI_gradshapeJ[f][j][i]);
            double aij = -0.5*D_avg*gij;

            MatSetValue(Aref,ir*G+gr,jr*G+gr,aij,ADD_VALUES);
          }//for j
        }//for i


//        // - Di^+ bj^-
//        for (int imap=0; imap<adj_fe_view->dofs; imap++)
//        {
//          int irmap = adj_ip_view->MapDof(imap);
//
//          for (int fj=0; fj<num_face_dofs; fj++)
//          {
//            int jmap  = MapCellDof(adj_cell,cell->faces[f]->v_indices[fj]);
//            int j     = MapCellDof(cell,cell->faces[f]->v_indices[fj]);
//            int jr    = cell_ip_view->MapDof(j);
//
//            double gij =
//              n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][jmap][imap]);
//            double aij = -0.5*adj_D[jmap]*gij;
//
//            MatSetValue(Aref,irmap*G+gr,jr*G+gr,aij,ADD_VALUES);
//          }//for j
//        }//for i

        //+ Di^- bj^+
        for (int fj=0; fj<num_face_dofs; fj++)
        {
          int j     = MapCellDof(cell, cell->faces[f].vertex_ids[fj]);
          int jmap  = MapCellDof(adj_cell, cell->faces[f].vertex_ids[fj]);
          int jrmap = adj_ip_view->MapDof(jmap);

          for (int i=0; i<fe_view->dofs; i++)
          {
            int ir = cell_ip_view->MapDof(i);

            double gij =
              n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][j][i]);
            double aij = 0.5*D_avg*gij;

            MatSetValue(Aref,ir*G+gr,jrmap*G+gr,aij,ADD_VALUES);
          }//for i
        }//for fj


        // - Dj^+ bi^-
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int imap  = MapCellDof(adj_cell, cell->faces[f].vertex_ids[fi]);
          int i     = MapCellDof(cell, cell->faces[f].vertex_ids[fi]);
          int ir    = cell_ip_view->MapDof(i);

          for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
          {
            int jrmap = adj_ip_view->MapDof(jmap);

            double gij =
              n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][imap][jmap]);
            double aij = -0.5*adj_D_avg*gij;

            MatSetValue(Aref,ir*G+gr,jrmap*G+gr,aij,ADD_VALUES);
          }//for j
        }//for i


      }//if not bndry
      else
      {
        int ir_boundary_index =
          abs(cell->faces[f].neighbor) - 1;
        int ir_boundary_type  = boundaries[ir_boundary_index]->type;

        if (ir_boundary_type == DIFFUSION_DIRICHLET)
        {
          //========================= Compute penalty coefficient
          double hm = HPerpendicular(cell, fe_view, f);

          //========================= Compute surface average D
          double D_avg = 0.0;
          double intS = 0.0;
          for (int fi=0; fi<num_face_dofs; fi++)
          {
            int i = fe_view->face_dof_mappings[f][fi];
            D_avg += D[i]*fe_view->IntS_shapeI[i][f];
            intS += fe_view->IntS_shapeI[i][f];
          }
          D_avg /= intS;

          double kappa = 1.0;
          if (cell->Type() == chi_mesh::CellType::SLAB)
            kappa = fmax(4.0*(D_avg/hm),0.25);
          if (cell->Type() == chi_mesh::CellType::POLYGON)
            kappa = fmax(4.0*(D_avg/hm),0.25);
          if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
            kappa = fmax(8.0*(D_avg/hm),0.25);

          //========================= Assembly penalty terms
          for (int fi=0; fi<num_face_dofs; fi++)
          {
            int i  = fe_view->face_dof_mappings[f][fi];
            int ir = cell_ip_view->MapDof(i);

            for (int fj=0; fj<num_face_dofs; fj++)
            {
              int j  = fe_view->face_dof_mappings[f][fj];
              int jr = cell_ip_view->MapDof(j);

              double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

              MatSetValue(Aref,ir*G+gr,jr*G+gr, aij,ADD_VALUES);
            }//for fj

          }//for fi

          // -Di^- bj^- and
          // -Dj^- bi^-
          for (int i=0; i<fe_view->dofs; i++)
          {
            int ir = cell_ip_view->MapDof(i);

            for (int j=0; j<fe_view->dofs; j++)
            {
              int jr = cell_ip_view->MapDof(j);

              double gij =
                n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][i][j] +
                      fe_view->IntS_shapeI_gradshapeJ[f][j][i]);
              double aij = -0.5*D_avg*gij;

              MatSetValue(Aref,ir*G+gr,jr*G+gr,aij,ADD_VALUES);
            }//for j
          }//for i
        }//Dirichlet
        else if (ir_boundary_type == DIFFUSION_ROBIN)
        {
          auto robin_bndry =
            (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];

          for (int fi=0; fi<num_face_dofs; fi++)
          {
            int i  = fe_view->face_dof_mappings[f][fi];
            int ir =  cell_ip_view->MapDof(i);

            for (int fj=0; fj<num_face_dofs; fj++)
            {
              int j  = fe_view->face_dof_mappings[f][fj];
              int jr =  cell_ip_view->MapDof(j);

              double aij = robin_bndry->a*fe_view->IntS_shapeI_shapeJ[f][i][j];
              aij /= robin_bndry->b;

              MatSetValue(Aref,ir*G+gr,jr*G+gr, aij,ADD_VALUES);
            }//for fj

            double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
            aii /= robin_bndry->b;

            MatSetValue(Aref,ir*G+gr,ir*G+gr, aii,ADD_VALUES);
          }//for fi
        }//robin
      }
    }//for f
  }//for gr
}

//###################################################################
/**Assembles b PWLD for polygon cells.*/
void chi_diffusion::Solver::PWLD_Assemble_b_GAGG(
                                               int cell_glob_index,
                                               chi_mesh::Cell *cell,
                                               DiffusionIPCellView* cell_ip_view)
{
  auto fe_view = (CellPWLFEView*)pwl_sdm->MapFeViewL(cell->local_id);

  for (int gr=0; gr<G; gr++)
  {
    //====================================== Process material id
    int mat_id = cell->material_id;

    std::vector<double> D(fe_view->dofs,1.0);
    std::vector<double> q(fe_view->dofs,1.0);
    std::vector<double> siga(fe_view->dofs,1.0);

    GetMaterialProperties(mat_id,cell,fe_view->dofs,D,q,siga,gi+gr);

    //========================================= Loop over DOFs
    for (int i=0; i<fe_view->dofs; i++)
    {
      int ir = cell_ip_view->MapDof(i);
      int ig = cell->vertex_ids[i];
      double rhsvalue =0.0;

      int ir_boundary_type;
      //if (!ApplyDirichletI(ir*G+gr,&ir_boundary_type,ig,true))
      {
        //====================== Develop matrix entry
        for (int j=0; j<fe_view->dofs; j++)
        {
          rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
        }//for j

        //====================== Apply RHS entry
        VecSetValue(bref,ir*G+gr,rhsvalue,ADD_VALUES);
      }//if ir not dirichlet

    }//for i
  }//for gr

}
