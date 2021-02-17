#include "diffusion_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/CellViews/pwl_cellbase.h"

#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_reflecting.h"
#include "Modules/DiffusionSolver/Boundaries/chi_diffusion_bndry_robin.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Assembles PWLC matrix for polygon cells.*/
void chi_diffusion::Solver::PWLD_Assemble_A_and_b(int cell_glob_index,
                                                  chi_mesh::Cell *cell,
                                                  DiffusionIPCellView* cell_ip_view,
                                                  int component,
                                                  int component_block_offset)
{
  auto pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWL>(this->discretization);
  auto fe_view = (CellPWLFEValues*)pwl_sdm->MapFeViewL(cell->local_id);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,0.0);

  GetMaterialProperties(mat_id, cell, fe_view->dofs, D, q, siga, component);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);
    double rhsvalue =0.0;

    //====================== Develop matrix entry
    for (int j=0; j<fe_view->dofs; j++)
    {
      int jr = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);

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
    auto& face = cell->faces[f];

    //================================== Get face normal
    chi_mesh::Vector3 n  = face.normal;

    int num_face_dofs = face.vertex_ids.size();

    if (face.has_neighbor)
    {
      chi_mesh::Cell*           adj_cell    = nullptr;
      CellPWLFEValues*               adj_fe_view = nullptr;
      int                              fmap = -1;

      //========================= Get adj cell information
      if (cell->faces[f].IsNeighborLocal(*grid))  //Local
      {
        adj_cell      = &grid->local_cells[face.GetNeighborLocalID(*grid)];
        adj_fe_view   = (CellPWLFEValues*)pwl_sdm->MapFeViewL(adj_cell->local_id);
      }//local
      else //Non-local
      {
        adj_cell    = pwl_sdm->MapNeighborCell(face.neighbor_id);
        adj_fe_view = pwl_sdm->MapNeighborCellFeView(face.neighbor_id);
      }//non-local
      //========================= Check valid information
      if (adj_cell == nullptr || adj_fe_view == nullptr)
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
                            adj_cell,
                            adj_fe_view->dofs,
                            adj_D,
                            adj_Q,
                            adj_sigma,
                            component);

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
        int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

        for (int fj=0; fj<num_face_dofs; fj++)
        {
          int j     = fe_view->face_dof_mappings[f][fj];
          int jr    = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);
          int jmap  = MapCellDof(adj_cell,cell->faces[f].vertex_ids[fj]);
          int jrmap = pwl_sdm->MapDOF(*adj_cell, jmap, unknown_manager, 0, component);

          double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

          MatSetValue(Aref,ir    ,jr   , aij,ADD_VALUES);
          MatSetValue(Aref,ir    ,jrmap,-aij,ADD_VALUES);
        }//for fj

      }//for fi


      //TODO: RESTORE POINT


      //========================= Assemble gradient terms
      // For the following comments we use the notation:
      // Dk = 0.5* n dot nabla bk

      // -Di^- bj^- and
      // -Dj^- bi^-
      for (int i=0; i<fe_view->dofs; i++)
      {
        int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

        for (int j=0; j<fe_view->dofs; j++)
        {
          int jr = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);

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
        int jrmap = pwl_sdm->MapDOF(*adj_cell, jmap, unknown_manager, 0, component);

        for (int i=0; i<fe_view->dofs; i++)
        {
          int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

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
        int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

        for (int jmap=0; jmap<adj_fe_view->dofs; jmap++)
        {
          int jrmap = pwl_sdm->MapDOF(*adj_cell, jmap, unknown_manager, 0, component);

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

//TODO: END RESTORE POINT

    }//if not bndry
    else
    {
      int ir_boundary_index = cell->faces[f].neighbor_id;
      int ir_boundary_type  = boundaries[ir_boundary_index]->type;

      if (ir_boundary_type == DIFFUSION_DIRICHLET)
      {
        auto dc_boundary =
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
          int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f][fj];
            int jr = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);

            double aij = kappa*fe_view->IntS_shapeI_shapeJ[f][i][j];

            MatSetValue(Aref,ir    ,jr, aij,ADD_VALUES);
            VecSetValue(bref,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
          }//for fj
        }//for fi

        // -Di^- bj^- and
        // -Dj^- bi^-
        for (int i=0; i<fe_view->dofs; i++)
        {
          int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

          for (int j=0; j<fe_view->dofs; j++)
          {
            int jr = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);

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
          int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_view->face_dof_mappings[f][fj];
            int jr = pwl_sdm->MapDOF(*cell, j, unknown_manager, 0, component);

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
                                            int component,
                                            int component_block_offset)
{
  auto pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWL>(this->discretization);
  auto fe_view = (CellPWLFEValues*)pwl_sdm->MapFeViewL(cell->local_id);

  //====================================== Process material id
  int mat_id = cell->material_id;

  std::vector<double> D(fe_view->dofs,1.0);
  std::vector<double> q(fe_view->dofs,1.0);
  std::vector<double> siga(fe_view->dofs,1.0);

  GetMaterialProperties(mat_id, cell, fe_view->dofs, D, q, siga, component);

  //========================================= Loop over DOFs
  for (int i=0; i<fe_view->dofs; i++)
  {
    int ir = pwl_sdm->MapDOF(*cell, i, unknown_manager, 0, component);

    //====================== Develop rhs entry
    double rhsvalue =0.0;
    for (int j=0; j<fe_view->dofs; j++)
      rhsvalue += q[j]*fe_view->IntV_shapeI_shapeJ[i][j];

    //====================== Apply RHS entry
    VecSetValue(bref,ir,rhsvalue,ADD_VALUES);

  }//for i


}
