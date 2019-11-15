#include "diffusion_solver.h"

#include <PiecewiseLinear/CellViews/pwl_cellbase.h>

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"
#include "../Boundaries/chi_diffusion_bndry_robin.h"

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::PWLD_Assemble_A_and_b(int cell_glob_index,
                                                  chi_mesh::Cell *cell,
                                                  DiffusionIPCellView* cell_ip_view,
                                                  int group)
{
  auto fe_view = (CellFEView*)pwl_discr->MapFeView(cell_glob_index);

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
    double rhsvalue =0.0;

    //====================== Develop matrix entry
    for (int j=0; j<fe_view->dofs; j++)
    {
      int jr =  cell_ip_view->MapDof(j);
      double jr_mat_entry =
        D[j]*fe_view->IntV_gradShapeI_gradShapeJ[i][j];

      jr_mat_entry +=
        siga[j]*fe_view->IntV_shapeI_shapeJ[i][j];

      MatSetValue(Aref,ir,jr,jr_mat_entry,ADD_VALUES);

      rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
    }//for j

    //====================== Apply RHS entry
    VecSetValue(bref,ir,rhsvalue,ADD_VALUES);

  }//for i


  //========================================= Loop over faces
  int num_faces = cell->faces.size();
  for (int f=0; f<num_faces; f++)
  {
    int neighbor = cell->faces[f].neighbor;

    //================================== Get face normal
    chi_mesh::Vector n  = cell->faces[f].normal;

    int num_face_dofs = cell->faces[f].vertex_ids.size();

    if (neighbor >=0)
    {
      try{
      chi_mesh::Cell*           adj_cell    = nullptr;
      CellFEView*               adj_fe_view = nullptr;
      DiffusionIPCellView*      adj_ip_view = nullptr;
      int                              fmap = -1;

      //========================= Get adj cell information
      if (grid->IsCellLocal(neighbor))  //Local
      {
        int adj_cell_local_index = grid->glob_cell_local_indices[neighbor];
        adj_ip_view   = ip_cell_views[adj_cell_local_index];
        adj_cell      = (chi_mesh::Cell*)grid->cells[neighbor];
        adj_fe_view   = (CellFEView*)pwl_discr->MapFeView(neighbor);
      }//local
      else //Non-local
      {
        int locI = grid->cells[neighbor]->partition_id;
        adj_ip_view = GetBorderIPView(locI,neighbor);
        adj_cell    = (chi_mesh::Cell*)GetBorderCell(locI,neighbor);
        adj_fe_view = (CellFEView*)GetBorderFEView(locI,neighbor);
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
      fmap = MapCellFace(cell,adj_cell,f);

      //========================= Compute penalty coefficient
      double hp = HPerpendicular(adj_cell, adj_fe_view, fmap);
      double hm = HPerpendicular(cell, fe_view, f);

      std::vector<double> adj_D,adj_Q,adj_sigma;

      GetMaterialProperties(adj_cell->material_id,
                            neighbor,
                            adj_fe_view->dofs,
                            adj_D,
                            adj_Q,
                            adj_sigma,
                            group);

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
        int imap = MapCellDof(adj_cell,cell->vertex_ids[i]);
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
          int j     = fe_view->face_dof_mappings[f][fj];
          int jr    = cell_ip_view->MapDof(j);
          int jmap  = MapCellDof(adj_cell,cell->faces[f].vertex_ids[fj]);
          int jrmap = adj_ip_view->MapDof(jmap);

          double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

          MatSetValue(Aref,ir    ,jr   , aij,ADD_VALUES);
          MatSetValue(Aref,ir    ,jrmap,-aij,ADD_VALUES);
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

          MatSetValue(Aref,ir,jr,aij,ADD_VALUES);
        }//for j
      }//for i


//      // - Di^+ bj^-
//      for (int imap=0; imap<adj_fe_view->dofs; imap++)
//      {
//        int irmap = adj_ip_view->MapDof(imap);
//
//        for (int fj=0; fj<num_face_dofs; fj++)
//        {
//          int jmap  = MapCellDof(adj_cell,polyh_cell->faces[f]->v_indices[fj]);
//          int j     = MapCellDof(polyh_cell,polyh_cell->faces[f]->v_indices[fj]);
//          int jr    = cell_ip_view->MapDof(j);
//
//          double gij =
//            n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][jmap][imap]);
//          double aij = -0.5*adj_D_avg*gij;
//
//          MatSetValue(Aref,irmap,jr,aij,ADD_VALUES);
//        }//for j
//      }//for i

      //+ Di^- bj^+
      for (int fj=0; fj<num_face_dofs; fj++)
      {
        int j     = MapCellDof(cell,cell->faces[f].vertex_ids[fj]);
        int jmap  = MapCellDof(adj_cell,cell->faces[f].vertex_ids[fj]);
        int jrmap = adj_ip_view->MapDof(jmap);

        for (int i=0; i<fe_view->dofs; i++)
        {
          int ir = cell_ip_view->MapDof(i);

          double gij =
            n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][j][i]);
          double aij = 0.5*D_avg*gij;

          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
        }//for i
      }//for fj



      // - Dj^+ bi^-
      for (int fi=0; fi<num_face_dofs; fi++)
      {
        int imap  = MapCellDof(adj_cell,cell->faces[f].vertex_ids[fi]);
        int i     = MapCellDof(cell,cell->faces[f].vertex_ids[fi]);
        int ir    = cell_ip_view->MapDof(i);

        for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
        {
          int jrmap = adj_ip_view->MapDof(jmap);

          double gij =
            n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][imap][jmap]);
          double aij = -0.5*adj_D_avg*gij;

          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
        }//for j
      }//for i

//      for (int fi=0; fi<num_face_dofs; fi++)
//      {
//        int i    = MapCellDof(polyh_cell,polyh_cell->faces[f]->v_indices[fi]);
//        int ir   = cell_ip_view->MapDof(i);
//        int imap = MapCellDof(adj_cell,polyh_cell->faces[f]->v_indices[fi]);
//
//        for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
//        {
//          int jrmap = adj_ip_view->MapDof(jmap);
//
//          double gij =
//            n.Dot(adj_fe_view->IntS_shapeI_gradshapeJ[fmap][jmap][imap]);
//          double aij = -0.5*adj_D_avg*gij;
//
//          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
//        }
//      }
//
//      for (int i=0; i<fe_view->dofs; i++)
//      {
//        int ir   = cell_ip_view->MapDof(i);
//
//        for (int fj=0; fj<num_face_dofs; fj++)
//        {
//          int j     = MapCellDof(polyh_cell,polyh_cell->faces[f]->v_indices[fj]);
//          int jmap  = MapCellDof(adj_cell,polyh_cell->faces[f]->v_indices[fj]);
//          int jrmap = adj_ip_view->MapDof(jmap);
//
//          double gij =
//            n.Dot(fe_view->IntS_shapeI_gradshapeJ[f][i][j]);
//          double aij = 0.5*D_avg*gij;
//
//          MatSetValue(Aref,ir,jrmap,aij,ADD_VALUES);
//        }
//      }

      }catch (std::out_of_range& o) {chi_log.Log(LOG_ALLERROR) << "OOR error.";exit(EXIT_FAILURE);}
    }//if not bndry
    else
    {
      int ir_boundary_index = abs(cell->faces[f].neighbor)-1;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_DIRICHLET)
      {
        chi_diffusion::BoundaryDirichlet* dc_boundary =
          (chi_diffusion::BoundaryDirichlet*)boundaries[ir_boundary_index];

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

            MatSetValue(Aref,ir    ,jr, aij,ADD_VALUES);
            VecSetValue(bref,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
          }//for fj

//          double rhs = kappa*fe_view->IntS_shapeI[i][f];
//          VecSetValue(bref,ir,rhs*0.1,ADD_VALUES);
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

            MatSetValue(Aref,ir,jr,aij,ADD_VALUES);
            VecSetValue(bref,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
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

            MatSetValue(Aref,ir ,jr, aij,ADD_VALUES);
          }//for fj

          double aii = robin_bndry->f*fe_view->IntS_shapeI[i][f];
          aii /= robin_bndry->b;

          MatSetValue(Aref,ir ,ir, aii,ADD_VALUES);
        }//for fi
      }//robin
    }
  }//for f
}

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::PWLD_Assemble_b(int cell_glob_index,
                                            chi_mesh::Cell *cell,
                                            DiffusionIPCellView* cell_ip_view,
                                            int group)
{
  auto fe_view = (CellFEView*)pwl_discr->MapFeView(cell_glob_index);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,1.0);

  GetMaterialProperties(mat_id,cell_glob_index,fe_view->dofs,D,q,siga,group);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = cell_ip_view->MapDof(i);
    int ig = cell->vertex_ids[i];
    double rhsvalue =0.0;

    int ir_boundary_type;
    //if (!ApplyDirichletI(ir,&ir_boundary_type,ig,true))
    {
      //====================== Develop matrix entry
      for (int j=0; j<fe_view->dofs; j++)
      {
        rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];
      }//for j

      //====================== Apply RHS entry
      VecSetValue(bref,ir,rhsvalue,ADD_VALUES);
    }//if ir not dirichlet

  }//for i


}