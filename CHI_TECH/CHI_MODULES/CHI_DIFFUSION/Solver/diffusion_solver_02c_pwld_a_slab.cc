#include "diffusion_solver.h"

#include <CHI_DISCRETIZATION_PWL/CellViews/pwl_slab.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

//###################################################################
/**Assembles A and b PWLD for slab cells.*/
void chi_diffusion::Solver::PWLD_Ab_Slab(int cell_glob_index,
                                            chi_mesh::Cell *cell,
                                            DIFFUSION_IP_VIEW* cell_ip_view,
                                            int group)
{
  chi_mesh::CellSlab* slab_cell =
    (chi_mesh::CellSlab*)(cell);
  SlabFEView* fe_view =
    (SlabFEView*)pwl_discr->MapFeView(cell_glob_index);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,group);


  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = cell_ip_view->MapDof(i);
    int ig = slab_cell->v_indices[i];
    double rhsvalue =0.0;

    int ir_boundary_type;
    //if (!ApplyDirichletI(ir,&ir_boundary_type,ig))
    {
      //====================== Develop matrix entry
      for (int j=0; j<fe_view->dofs; j++)
      {
        int jr =  cell_ip_view->MapDof(j);
        int jg = slab_cell->v_indices[j];
        double jr_mat_entry =
          D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

        jr_mat_entry +=
          siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

        int jr_boundary_type;
        //if (!ApplyDirichletJ(jr,ir,jr_mat_entry,&jr_boundary_type,jg))
        {
          MatSetValue(Aref,ir,jr,jr_mat_entry,ADD_VALUES);
        }

        rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
      }//for j

      //====================== Develop RHS entry
      VecSetValue(bref,ir,rhsvalue,ADD_VALUES);
    }//if ir not dirichlet

  }//for i

  //========================================= Loop over faces
  int num_faces = 2;
  for (int f=0; f<num_faces; f++)
  {
    int neighbor = slab_cell->edges[f];

    //================================== Get face normal
    chi_mesh::Vertex v0 = *grid->nodes[slab_cell->v_indices[0]];
    chi_mesh::Vertex v1 = *grid->nodes[slab_cell->v_indices[1]];
    chi_mesh::Vector n = slab_cell->face_normals[f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERFACE TERMS
    if (neighbor >=0)
    {
      chi_mesh::CellSlab*    adj_cell = nullptr;
      SlabFEView*         adj_fe_view = nullptr;
      DIFFUSION_IP_VIEW*  adj_ip_view = nullptr;

      //========================= Get adj cell information
      if (grid->IsCellLocal(neighbor))  //Local
      {
        int adj_cell_local_index = grid->glob_cell_local_indices[neighbor];
        adj_ip_view   = ip_cell_views[adj_cell_local_index];
        adj_cell      = (chi_mesh::CellSlab*)grid->cells[neighbor];
        adj_fe_view   = (SlabFEView*)pwl_discr->MapFeView(neighbor);
      }//local
      else //Non-local
      {
        int locI = grid->cells[neighbor]->partition_id;
        adj_ip_view = GetBorderIPView(locI,neighbor);
        adj_cell    = (chi_mesh::CellSlab*)GetBorderCell(locI,neighbor);
        adj_fe_view = (SlabFEView*)GetBorderFEView(locI,neighbor);
      }//non-local

      //========================= Check valid information
      if (adj_cell == nullptr || adj_fe_view == nullptr ||
          adj_ip_view == nullptr)
      {
        chi_log.Log(LOG_ALL)
          << "Error in MIP cell information.";
        exit(EXIT_FAILURE);
      }

      //========================= Compute penalty coefficient
      double hp = adj_fe_view->h/2.0;
      double hm = fe_view->h/2.0;
//      double adj_D,adj_Q,adj_sigma;
//      GetMaterialProperties(adj_cell->material_id,adj_D,adj_Q,adj_sigma);

      std::vector<double> adj_D,adj_Q,adj_sigma;

      GetMaterialProperties(adj_cell->material_id,
                            neighbor,
                            adj_fe_view->dofs,
                            adj_D,
                            adj_Q,
                            adj_sigma,
                            group);

      double kappa = fmax(2.0*(adj_D[f]/hp + D[f]/hm),0.25);

      //========================= Assembly penalty terms
      for (int fi=f; fi<(f+1); fi++)
      {
        int i  = fi;
        int ir = cell_ip_view->MapDof(i);

        //Mapping face index to adj-cell
        int imap = MapCellDof(adj_cell,slab_cell->v_indices[fi]);
        int irstar = adj_ip_view->MapDof(imap);

        for (int fj=f; fj<(f+1); fj++)
        {
          int j  = fj;
          int jr = cell_ip_view->MapDof(j);
          int jmap  = MapCellDof(adj_cell,slab_cell->v_indices[fj]);
          int jrmap = adj_ip_view->MapDof(jmap);

          double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];


          MatSetValue(Aref,ir    ,jr   , aij,ADD_VALUES);
          MatSetValue(Aref,ir    ,jrmap,-aij,ADD_VALUES);
        }//for fj

      }//for fi

      //========================= Assemble gradient terms
      //Map the current face to the adjacent cell's face number
      int fmap = -1;
      for (int af=0; af<2; af++)
        if (slab_cell->v_indices[f] == adj_cell->v_indices[af])
          fmap = af;

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
          double aij = -0.5*D[f]*gij;

          MatSetValue(Aref,ir,jr,aij,ADD_VALUES);
        }//for j
      }//for i

      // - Di^+ bj^-
//      for (int imap=0; imap<adj_fe_view->dofs; imap++)
//      {
//        int irmap = adj_ip_view->MapDof(imap);
//
//        for (int fj=f; fj<(f+1); fj++)
//        {
//          int jmap  = MapCellDof(adj_cell,slab_cell->v_indices[fj]);
//          int j     = MapCellDof(slab_cell,slab_cell->v_indices[fj]);
//          int jr    = cell_ip_view->MapDof(j);
//
//          double gij =
//            n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][jmap][imap]);
//          double aij = -0.5*adj_D[f]*gij;
//
//          MatSetValue(Aref,irmap,jr,aij,ADD_VALUES);
//        }//for j
//      }//for i

      //+ Di^- bj^+
      for (int fj=0; fj<2; fj++)
      {
        int j     = MapCellDof(slab_cell,slab_cell->v_indices[fj]);
        int jmap  = MapCellDof(adj_cell,slab_cell->v_indices[fj]);
        int jrmap = adj_ip_view->MapDof(jmap);

        for (int i=0; i<fe_view->dofs; i++)
        {
          int ir = cell_ip_view->MapDof(i);

          double gij =
            n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][j][i]);
          double aij = 0.5*D[f]*gij;

          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
        }//for i
      }//for fj

      // - Dj^+ bi^-
      for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
      {
        int jrmap = adj_ip_view->MapDof(jmap);

        for (int fi=f; fi<(f+1); fi++)
        {
          int imap  = MapCellDof(adj_cell,slab_cell->v_indices[fi]);
          int i     = MapCellDof(slab_cell,slab_cell->v_indices[fi]);
          int ir    = cell_ip_view->MapDof(i);

          double gij =
            n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][imap][jmap]);
          double aij = -0.5*adj_D[f]*gij;

          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
        }//for j
      }//for i

    }//if not bndry
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOUNDARY TERMS
    else
    {
      int ig = slab_cell->v_indices[f];
      int ir_boundary_index = abs(nodal_boundary_numbers[ig])-1;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_DIRICHLET)
      {
        chi_diffusion::BoundaryDirichlet* dc_boundary =
          (chi_diffusion::BoundaryDirichlet*)boundaries[ir_boundary_index];

        //========================= Compute penalty coefficient
        double hm = fe_view->h/2.0;
        double kappa = fmax(4.0*D[f]/hm,0.25);

        //========================= Assembly penalty terms
        for (int fi=f; fi<(f+1); fi++)
        {
          int i  = fi;
          int ir = cell_ip_view->MapDof(i);

          for (int fj=f; fj<(f+1); fj++)
          {
            int j  = fj;
            int jr = cell_ip_view->MapDof(j);

            double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

            MatSetValue(Aref,ir    ,jr, aij,ADD_VALUES);
            VecSetValue(bref,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
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
            double aij = -0.5*D[f]*gij;

            MatSetValue(Aref,ir,jr,aij,ADD_VALUES);
            VecSetValue(bref,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
          }//for j
        }//for i

      }//if dirichlet
      else if (ir_boundary_type == DIFFUSION_ROBIN)
      {
        auto robin_bndry =
          (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];

        for (int fi=f; fi<(f+1); fi++)
        {
          int i  = fi;
          int ir = cell_ip_view->MapDof(i);

          for (int fj=f; fj<(f+1); fj++)
          {
            int j  = fj;
            int jr = cell_ip_view->MapDof(j);

            double aij = robin_bndry->a*fe_view->IntS_shapeI_shapeJ[f][i][j];
            aij /= robin_bndry->b;

            MatSetValue(Aref,ir ,jr, aij,ADD_VALUES);
          }//for fj

          double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
          aii /= robin_bndry->b;

          MatSetValue(Aref,ir ,ir, aii,ADD_VALUES);
        }//for fi
      }//robin
    }//if bndry

  }//for f
}


//###################################################################
/**Assembles b PWLD for slab cells.*/
void chi_diffusion::Solver::PWLD_b_Slab(int cell_glob_index,
                                         chi_mesh::Cell *cell,
                                         DIFFUSION_IP_VIEW* cell_ip_view,
                                         int group)
{
  chi_mesh::CellSlab* slab_cell =
    (chi_mesh::CellSlab*)(cell);
  SlabFEView* fe_view =
    (SlabFEView*)pwl_discr->MapFeView(cell_glob_index);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,group);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = cell_ip_view->MapDof(i);
    int ig = slab_cell->v_indices[i];
    double rhsvalue =0.0;

    int ir_boundary_type;
    //if (!ApplyDirichletI(ir,&ir_boundary_type,ig))
    {
      //====================== Develop matrix entry
      for (int j=0; j<fe_view->dofs; j++)
      {
        rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
      }//for j

      //====================== Develop RHS entry
      VecSetValue(bref,ir,rhsvalue,ADD_VALUES);
    }//if ir not dirichlet

  }//for i
}