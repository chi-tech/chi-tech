#include "diffusion_solver.h"
//
//#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"
//
//#include "chi_runtime.h"
//
////###################################################################
///**Assembles PWLC matrix for polygon cells.*/
//void chi_diffusion::Solver::PWLD_Assemble_A_and_b_GAGG(const chi_mesh::Cell& cell)
//{
//  auto pwl_sdm = std::static_pointer_cast<chi_math::SpatialDiscretization_PWLD>(this->discretization);
////  auto fe_view = (CellPWLFEValues*)pwl_sdm->MapFeViewL(cell.local_id);
//  const auto& fe_intgrl_values = pwl_sdm->GetUnitIntegrals(cell);
//
//  size_t num_nodes = fe_intgrl_values.NumNodes();
//
//  for (int gr=0; gr<G; gr++)
//  {
//    //====================================== Process material id
//    int mat_id = cell.material_id;
//
//    std::vector<double> D(num_nodes, 1.0);
//    std::vector<double> q(num_nodes, 1.0);
//    std::vector<double> siga(num_nodes, 1.0);
//
//    GetMaterialProperties(cell, num_nodes, D, q, siga, gi + gr);
//
//    //========================================= Loop over DOFs
//    for (int i=0; i<num_nodes; i++)
//    {
//      int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//      double rhsvalue =0.0;
//
//      //====================== Develop matrix entry
//      for (int j=0; j<num_nodes; j++)
//      {
//        int jr = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//
//        double jr_mat_entry =
//          D[j]* fe_intgrl_values.IntV_gradShapeI_gradShapeJ(i, j);
//
//        jr_mat_entry +=
//          siga[j]* fe_intgrl_values.IntV_shapeI_shapeJ(i, j);
//
//        MatSetValue(A,ir,jr,jr_mat_entry,ADD_VALUES);
//
//        rhsvalue += q[j]* fe_intgrl_values.IntV_shapeI_shapeJ(i, j);
//      }//for j
//
//      //====================== Apply RHS entry
//      VecSetValue(b,ir,rhsvalue,ADD_VALUES);
//
//    }//for i
//
//    //========================================= Loop over faces
//    int num_faces = cell.faces.size();
//    for (unsigned int f=0; f<num_faces; f++)
//    {
//      auto& face = cell.faces[f];
//
//      //================================== Get face normal
//      chi_mesh::Vector3 n  = face.normal;
//
//      int num_face_dofs = face.vertex_ids.size();
//
//      if (face.has_neighbor)
//      {
//        const auto& adj_cell = grid->cells[face.neighbor_id];
//        const auto& adj_fe_intgrl_values = pwl_sdm->GetUnitIntegrals(adj_cell);
//
//        //========================= Get the current map to the adj cell's face
//        unsigned int fmap = MapCellFace(cell,adj_cell,f);
//
//        //========================= Compute penalty coefficient
//        double hp = HPerpendicular(adj_cell, adj_fe_intgrl_values, fmap);
//        double hm = HPerpendicular(cell, fe_intgrl_values, f);
//
//        std::vector<double> adj_D,adj_Q,adj_sigma;
//
//        GetMaterialProperties(adj_cell,
//                              adj_fe_intgrl_values.NumNodes(),
//                              adj_D,
//                              adj_Q,
//                              adj_sigma,
//                              gi+gr);
//
//        //========================= Compute surface average D
//        double D_avg = 0.0;
//        double intS = 0.0;
//        for (int fi=0; fi<num_face_dofs; fi++)
//        {
//          int i = fe_intgrl_values.FaceDofMapping(f,fi);
//          D_avg += D[i]* fe_intgrl_values.IntS_shapeI(f, i);
//          intS += fe_intgrl_values.IntS_shapeI(f, i);
//        }
//        D_avg /= intS;
//
//        //========================= Compute surface average D_adj
//        double adj_D_avg = 0.0;
//        double adj_intS = 0.0;
//        for (int fi=0; fi<num_face_dofs; fi++)
//        {
//          int i    = fe_intgrl_values.FaceDofMapping(f,fi);
//          int imap = MapCellLocalNodeIDFromGlobalID(adj_cell, cell.vertex_ids[i]);
//          adj_D_avg += adj_D[imap]* adj_fe_intgrl_values.IntS_shapeI(fmap, imap);
//          adj_intS += adj_fe_intgrl_values.IntS_shapeI(fmap, imap);
//        }
//        adj_D_avg /= adj_intS;
//
//        //========================= Compute kappa
//        double kappa = 1.0;
//        if (cell.Type() == chi_mesh::CellType::SLAB)
//          kappa = fmax(2.0*(adj_D_avg/hp + D_avg/hm),0.25);
//        if (cell.Type() == chi_mesh::CellType::POLYGON)
//          kappa = fmax(2.0*(adj_D_avg/hp + D_avg/hm),0.25);
//        if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
//          kappa = fmax(4.0*(adj_D_avg/hp + D_avg/hm),0.25);
//
//        //========================= Assembly penalty terms
//        for (int fi=0; fi<num_face_dofs; fi++)
//        {
//          int i  = fe_intgrl_values.FaceDofMapping(f,fi);
//          int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//
//          for (int fj=0; fj<num_face_dofs; fj++)
//          {
//            int j     = fe_intgrl_values.FaceDofMapping(f,fj);
//            int jr    = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//            int jmap  = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids[fj]);
//            int jrmap = pwl_sdm->MapDOF(adj_cell, jmap, unknown_manager, 0, gr);
//
//            double aij = kappa* fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
//
//            MatSetValue(A,ir    ,jr   , aij,ADD_VALUES);
//            MatSetValue(A,ir    ,jrmap,-aij,ADD_VALUES);
//          }//for fj
//
//        }//for fi
//
//        //========================= Assemble gradient terms
//        // For the following comments we use the notation:
//        // Dk = 0.5* n dot nabla bk
//
//        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
//        for (int i=0; i<fe_intgrl_values.NumNodes(); i++)
//        {
//          int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//
//          for (int fj=0; fj<num_face_dofs; fj++)
//          {
//            int j     = fe_intgrl_values.FaceDofMapping(f,fj);
//            int jr    = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//            int jmap  = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids[fj]);
//            int jrmap = pwl_sdm->MapDOF(adj_cell, jmap, unknown_manager, 0, gr);
//
//            double aij =
//              -0.5*D_avg*n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, j, i));
//
//            MatSetValue(A,ir,jr   , aij,ADD_VALUES);
//            MatSetValue(A,ir,jrmap,-aij,ADD_VALUES);
//          }//for fj
//        }//for i
//
//        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
//        for (int fi=0; fi<num_face_dofs; fi++)
//        {
//          int i     = fe_intgrl_values.FaceDofMapping(f,fi);
//          int ir    = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//          int imap  = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids[fi]);
//          int irmap = pwl_sdm->MapDOF(adj_cell, imap, unknown_manager, 0, gr);
//
//          for (int j=0; j<fe_intgrl_values.NumNodes(); j++)
//          {
//            int jr = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//
//            double aij =
//              -0.5*D_avg*n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, i, j));
//
//            MatSetValue(A,ir   ,jr, aij,ADD_VALUES);
//            MatSetValue(A,irmap,jr,-aij,ADD_VALUES);
//          }//for j
//        }//for fi
//
//
//      }//if not bndry
//      else
//      {
//        int ir_boundary_index = face.neighbor_id;
//        auto ir_boundary_type  = boundaries[ir_boundary_index]->type;
//
//        if (ir_boundary_type == BoundaryType::Dirichlet)
//        {
////          auto dc_boundary =
////            (chi_diffusion::BoundaryDirichlet*)boundaries[ir_boundary_index];
//
//          //========================= Compute penalty coefficient
//          double hm = HPerpendicular(cell, fe_intgrl_values, f);
//
//          //========================= Compute surface average D
//          double D_avg = 0.0;
//          double intS = 0.0;
//          for (int fi=0; fi<num_face_dofs; fi++)
//          {
//            int i = fe_intgrl_values.FaceDofMapping(f,fi);
//            D_avg += D[i]* fe_intgrl_values.IntS_shapeI(f, i);
//            intS += fe_intgrl_values.IntS_shapeI(f, i);
//          }
//          D_avg /= intS;
//
//          double kappa = 1.0;
//          if (cell.Type() == chi_mesh::CellType::SLAB)
//            kappa = fmax(4.0*(D_avg/hm),0.25);
//          if (cell.Type() == chi_mesh::CellType::POLYGON)
//            kappa = fmax(4.0*(D_avg/hm),0.25);
//          if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
//            kappa = fmax(8.0*(D_avg/hm),0.25);
//
//          //========================= Assembly penalty terms
//          for (int fi=0; fi<num_face_dofs; fi++)
//          {
//            int i  = fe_intgrl_values.FaceDofMapping(f,fi);
//            int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//
//            for (int fj=0; fj<num_face_dofs; fj++)
//            {
//              int j  = fe_intgrl_values.FaceDofMapping(f,fj);
//              int jr = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//
//              double aij = kappa* fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
//
//              MatSetValue(A,ir    ,jr, aij,ADD_VALUES);
////              VecSetValue(b,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
//            }//for fj
//          }//for fi
//
//          // -Di^- bj^- and
//          // -Dj^- bi^-
//          for (int i=0; i<fe_intgrl_values.NumNodes(); i++)
//          {
//            int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//
//            for (int j=0; j<fe_intgrl_values.NumNodes(); j++)
//            {
//              int jr = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//
//              double gij =
//                n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, i, j) +
//                      fe_intgrl_values.IntS_shapeI_gradshapeJ(f, j, i));
//              double aij = -0.5*D_avg*gij;
//
//              MatSetValue(A,ir,jr,aij,ADD_VALUES);
////              VecSetValue(b,ir,aij*dc_boundary->boundary_value,ADD_VALUES);
//            }//for j
//          }//for i
//        }//Dirichlet
//        else if (ir_boundary_type == BoundaryType::Robin)
//        {
//          auto robin_bndry =
//            (chi_diffusion::BoundaryRobin*)boundaries[ir_boundary_index];
//
//          for (int fi=0; fi<num_face_dofs; fi++)
//          {
//            int i  = fe_intgrl_values.FaceDofMapping(f,fi);
//            int ir = pwl_sdm->MapDOF(cell, i, unknown_manager, 0, gr);
//
//            for (int fj=0; fj<num_face_dofs; fj++)
//            {
//              int j  = fe_intgrl_values.FaceDofMapping(f,fj);
//              int jr = pwl_sdm->MapDOF(cell, j, unknown_manager, 0, gr);
//
//              double aij = robin_bndry->a* fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
//              aij /= robin_bndry->b;
//
//              MatSetValue(A,ir ,jr, aij,ADD_VALUES);
//            }//for fj
//
//            double aii = robin_bndry->f* fe_intgrl_values.IntS_shapeI(f, i);
//            aii /= robin_bndry->b;
//
//            MatSetValue(A,ir ,ir, aii,ADD_VALUES);
//          }//for fi
//        }//robin
//      }
//    }//for f
//  }//for gr
//}

////###################################################################
///**Assembles b PWLD for polygon cells.*/
//void chi_diffusion::Solver::PWLD_Assemble_b_GAGG(const chi_mesh::Cell& cell)
//{
//  auto pwl_sdm = std::static_pointer_cast<chi_math::SpatialDiscretization_PWLD>(this->discretization);
//  const auto& fe_intgrl_values = pwl_sdm->GetUnitIntegrals(cell);
//
//  size_t num_nodes = fe_intgrl_values.NumNodes();
//
//  for (int gr=0; gr<G; gr++)
//  {
//    //====================================== Process material
//    std::vector<double> D(num_nodes, 1.0);
//    std::vector<double> q(num_nodes, 1.0);
//    std::vector<double> siga(num_nodes, 1.0);
//
//    GetMaterialProperties(cell, num_nodes, D, q, siga, gi + gr);
//
//    //========================================= Loop over DOFs
//    for (int i=0; i<num_nodes; i++)
//    {
//      int ir = pwl_sdm->MapDOF(cell,i,unknown_manager,0,gr);
//
//      //====================== Develop rhs entry
//      double rhsvalue =0.0;
//      for (int j=0; j<num_nodes; j++)
//        rhsvalue += q[j]* fe_intgrl_values.IntV_shapeI_shapeJ(i, j);
//
//      //====================== Apply RHS entry
//      VecSetValue(b,ir,rhsvalue,ADD_VALUES);
//
//    }//for i
//  }//for gr
//
//}
