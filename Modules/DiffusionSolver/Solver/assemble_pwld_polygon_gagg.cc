#include "diffusion_solver.h"

#include <PiecewiseLinear/CellViews/pwl_polygon.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Assembles A and b PWLD for polygon cells.*/
void chi_diffusion::Solver::PWLD_Ab_Polygon_GAGG(int cell_glob_index,
                                            chi_mesh::Cell *cell,
                                            DiffusionIPCellView* cell_ip_view)
{
  chi_mesh::CellPolygon* poly_cell =
    (chi_mesh::CellPolygon*)(cell);
  PolygonFEView* fe_view =
    (PolygonFEView*)pwl_discr->MapFeView(cell_glob_index);

  for (int gr=0; gr<G; gr++)
  {
    //====================================== Process material id
    int mat_id = cell->material_id;

    std::vector<double> D(fe_view->dofs,1.0);
    std::vector<double> q(fe_view->dofs,1.0);
    std::vector<double> siga(fe_view->dofs,0.0);

    GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,gi+gr);

    //========================================= Loop over DOFs
    for (int i=0; i<fe_view->dofs; i++)
    {
      int ir = cell_ip_view->MapDof(i);
      int ig = poly_cell->v_indices[i];
      double rhsvalue =0.0;

      int ir_boundary_type;
      if (!ApplyDirichletI(ir,&ir_boundary_type,ig))
      {
        //====================== Develop matrix entry
        for (int j=0; j<fe_view->dofs; j++)
        {
          int jr =  cell_ip_view->MapDof(j);
          int jg = poly_cell->v_indices[j];
          double jr_mat_entry =
            D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

          jr_mat_entry +=
            siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

          int jr_boundary_type;
          if (!ApplyDirichletJ(jr,ir,jr_mat_entry,&jr_boundary_type,jg))
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
    int num_faces = poly_cell->edges.size();
    for (int f=0; f<num_faces; f++)
    {
      int neighbor = poly_cell->edges[f][2];

      //================================== Get face normal
      chi_mesh::Vertex v0 = *grid->nodes[poly_cell->edges[f][0]];
      chi_mesh::Vertex v1 = *grid->nodes[poly_cell->edges[f][1]];
      chi_mesh::Vector n  = (v1-v0).Cross(chi_mesh::Vector(0,0,1));
      n = n/n.Norm();

      int num_face_dofs = 2;

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERFACE TERMS
      if (neighbor >=0)
      {
        chi_mesh::CellPolygon* adj_cell = nullptr;
        PolygonFEView*         adj_fe_view = nullptr;
        DiffusionIPCellView*     adj_ip_view = nullptr;

        //========================= Get adj cell information
        if (grid->IsCellLocal(neighbor))  //Local
        {
          int adj_cell_local_index = grid->glob_cell_local_indices[neighbor];
          adj_ip_view   = ip_cell_views[adj_cell_local_index];
          adj_cell      = (chi_mesh::CellPolygon*)grid->cells[neighbor];
          adj_fe_view   = (PolygonFEView*)pwl_discr->MapFeView(neighbor);
        }//local
        else //Non-local
        {
          int locI = grid->cells[neighbor]->partition_id;
          adj_ip_view = GetBorderIPView(locI,neighbor);
          adj_cell    = (chi_mesh::CellPolygon*)GetBorderCell(locI,neighbor);
          adj_fe_view = (PolygonFEView*)GetBorderFEView(locI,neighbor);
        }//non-local

        //========================= Check valid information
        if (adj_cell == nullptr || adj_fe_view == nullptr ||
            adj_ip_view == nullptr)
        {
          chi_log.Log(LOG_ALL)
            << "Error in MIP cell information.";
          exit(EXIT_FAILURE);
        }

        //========================= Assemble gradient terms
        //Map the current face to the adjacent cell's face number
        int fmap = -1;
        for (int af=0; af<adj_cell->edges.size(); af++)
          if (poly_cell->edges[f][0] == adj_cell->edges[af][1])
            fmap = af;

        //========================= Compute penalty coefficient
        double area_m  = 0.0;
        for (int i=0; i<fe_view->dofs; i++)
          area_m += fe_view->IntV_shapeI[i];

        double area_p  = 0.0;
        for (int i=0; i<adj_fe_view->dofs; i++)
          area_p += adj_fe_view->IntV_shapeI[i];

        double hp = HPerpendicularPoly(4, area_p, (v1 - v0).Norm());
        double hm = HPerpendicularPoly(4, area_m, (v1 - v0).Norm());

        std::vector<double> adj_D,adj_Q,adj_sigma;

        GetMaterialProperties(adj_cell->material_id,
                              neighbor,
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
        for (int fi=0; fi<2; fi++)
        {
          int i    = fe_view->face_dof_mappings[f][fi];
          int imap = MapCellDof(adj_cell,poly_cell->v_indices[i]);
          adj_D_avg += adj_D[imap]*adj_fe_view->IntS_shapeI[imap][fmap];
          adj_intS += adj_fe_view->IntS_shapeI[imap][fmap];
        }
        adj_D_avg /= adj_intS;

        double kappa = fmax(2.0*(adj_D_avg/hp + D_avg/hm),0.25);

        //========================= Assembly penalty terms
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_view->face_dof_mappings[f][fi];
          int ir = cell_ip_view->MapDof(i);

          //Mapping face index to adj-cell
          int imap = MapCellDof(adj_cell,poly_cell->edges[f][fi]);
          int irstar = adj_ip_view->MapDof(imap);

          for (int fj=0; fj<2; fj++)
          {
            int j  = fe_view->face_dof_mappings[f][fj];
            int jr = cell_ip_view->MapDof(j);

            double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

            MatSetValue(Aref,ir*G+gr    ,jr*G+gr, aij,ADD_VALUES);
            MatSetValue(Aref,irstar*G+gr,jr*G+gr,-aij,ADD_VALUES);
          }//for fj

        }//for fi

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

        // - Di^+ bj^-
        for (int imap=0; imap<adj_fe_view->dofs; imap++)
        {
          int irmap = adj_ip_view->MapDof(imap);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int jmap  = MapCellDof(adj_cell,poly_cell->edges[f][fj]);
            int j     = MapCellDof(poly_cell,poly_cell->edges[f][fj]);
            int jr    = cell_ip_view->MapDof(j);

            double gij =
              n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][jmap][imap]);
            double aij = -0.5*adj_D_avg*gij;

            MatSetValue(Aref,irmap*G+gr,jr*G+gr,aij,ADD_VALUES);
          }//for j
        }//for i

        // - Dj^+ bi^-
        for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
        {
          int jrmap = adj_ip_view->MapDof(jmap);

          for (int fi=0; fi<num_face_dofs; fi++)
          {
            int imap  = MapCellDof(adj_cell,poly_cell->edges[f][fi]);
            int i     = MapCellDof(poly_cell,poly_cell->edges[f][fi]);
            int ir    = cell_ip_view->MapDof(i);

            double gij =
              n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][imap][jmap]);
            double aij = -0.5*adj_D_avg*gij;

            MatSetValue(Aref,ir*G+gr,jrmap*G+gr,aij,ADD_VALUES);
          }//for j
        }//for i

      }//if not bndry
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY TERMS
      else
      {
        int ir_boundary_index = abs(poly_cell->edges[f][2])-1;
        int ir_boundary_type  = boundaries[ir_boundary_index]->type;

        if (ir_boundary_type == DIFFUSION_DIRICHLET)
        {
          //========================= Compute penalty coefficient
          double area_m  = 0.0;
          for (int i=0; i<fe_view->dofs; i++)
            area_m += fe_view->IntV_shapeI[i];

          double hm = HPerpendicularPoly(4, area_m, (v1 - v0).Norm());

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

          double kappa = fmax(4.0*(D_avg/hm),0.25);

          //========================= Assembly penalty terms
          for (int fi=0; fi<num_face_dofs; fi++)
          {
            int i  = fe_view->face_dof_mappings[f][fi];
            int ir = cell_ip_view->MapDof(i);

            for (int fj=0; fj<2; fj++)
            {
              int j  = fe_view->face_dof_mappings[f][fj];
              int jr = cell_ip_view->MapDof(j);

              double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

              MatSetValue(Aref,ir*G+gr,jr*G+gr, aij,ADD_VALUES);
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
        }//if dirichlet
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
      }//if bndry
    }//for f
  }//for gr

}

//###################################################################
/**Assembles b PWLD for polygon cells.*/
void chi_diffusion::Solver::PWLD_b_Polygon_GAGG(
                                    int cell_glob_index,
                                    chi_mesh::Cell *cell,
                                    DiffusionIPCellView* cell_ip_view)
{
  chi_mesh::CellPolygon* poly_cell =
    (chi_mesh::CellPolygon*)(cell);
  PolygonFEView* fe_view =
    (PolygonFEView*)pwl_discr->MapFeView(cell_glob_index);

  for (int gr=0; gr<G; gr++)
  {
    //====================================== Process material id
    int mat_id = cell->material_id;

    std::vector<double> D(fe_view->dofs,1.0);
    std::vector<double> q(fe_view->dofs,1.0);
    std::vector<double> siga(fe_view->dofs,0.0);

    GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,gi+gr);

    //========================================= Loop over DOFs
    for (int i=0; i<fe_view->dofs; i++)
    {
      int ir = cell_ip_view->MapDof(i);
      int ig = poly_cell->v_indices[i];
      double rhsvalue =0.0;

      int ir_boundary_type;
      if (!ApplyDirichletI(ir,&ir_boundary_type,ig))
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