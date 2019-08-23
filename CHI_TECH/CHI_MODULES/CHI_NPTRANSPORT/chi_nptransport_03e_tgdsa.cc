#include "chi_nptransport.h"

#include <CHI_MESH/CHI_CELL/cell.h>
#include <CHI_MESH/CHI_CELL/cell_slab.h>
#include <CHI_MESH/CHI_CELL/cell_polygon.h>
#include <CHI_MESH/CHI_CELL/cell_polyhedron.h>

#include <CHI_MODULES/CHI_DIFFUSION/Solver/diffusion_solver.h>
#include <CHI_MODULES/CHI_DIFFUSION/Boundaries/chi_diffusion_bndry_dirichlet.h>
#include <CHI_PHYSICS/chi_physics.h>
#include <chi_log.h>

extern CHI_LOG chi_log;
extern CHI_PHYSICS chi_physics_handler;

//###################################################################
/**Initializes the Within-Group DSA solver. */
void CHI_NPTRANSPORT::InitTGDSA(NPT_GROUPSET *groupset)
{
  if (groupset->apply_tgdsa)
  {
    //================================= Initialize field function
    delta_phi_local.resize(local_dof_count,0.0);
    int g = 0;
    int m = 0;
    chi_physics::FieldFunction* deltaphi_ff =
      new chi_physics::FieldFunction;
    deltaphi_ff->text_name = std::string("Sum_Sigma_s_DeltaPhi_g") +
                             std::to_string(g) +
                             std::string("_m") + std::to_string(m);
    deltaphi_ff->grid = grid;
    deltaphi_ff->spatial_discretization = discretization;
    deltaphi_ff->id = chi_physics_handler.fieldfunc_stack.size();

    deltaphi_ff->type = FF_SDM_PWLD;
    deltaphi_ff->num_grps = 1;
    deltaphi_ff->num_moms = 1;
    deltaphi_ff->grp = g;
    deltaphi_ff->mom = m;
    deltaphi_ff->field_vector_local = &delta_phi_local;

    deltaphi_ff->local_cell_dof_array_address =
      &local_cell_dof_array_address;

    chi_physics_handler.fieldfunc_stack.push_back(deltaphi_ff);
    field_functions.push_back(deltaphi_ff);

    //================================= Set diffusion solver
    std::string solver_name = std::string("TGDSA");
    solver_name += std::string("[g=")
                 + std::to_string(groupset->groups.front()->id)
                 + std::string("-")
                 + std::to_string(groupset->groups.back()->id)
                 + std::string("]");
    chi_diffusion::Solver* dsolver = new chi_diffusion::Solver(solver_name);
    groupset->tgdsa_solver = dsolver;

    dsolver->regions.push_back(this->regions.back());
    dsolver->discretization = discretization;
    dsolver->fem_method = PWLD_MIP;
    dsolver->residual_tolerance = groupset->tgdsa_tol;
    dsolver->max_iters          = groupset->tgdsa_max_iters;
    dsolver->options_string     = groupset->tgdsa_string;
    if (groupset->apply_wgdsa)
      dsolver->material_mode = DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL;
    else
      dsolver->material_mode = DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART;
    dsolver->q_field = deltaphi_ff;

    //================================= Initialize boundaries
    if (not dsolver->common_items_initialized)
      dsolver->InitializeCommonItems();

    for (int b=0; b<dsolver->boundaries.size(); b++)
    {
//      chi_diffusion::BoundaryRobin* bound =
//        new chi_diffusion::BoundaryRobin(0.25,0.5,0.0);
      chi_diffusion::BoundaryDirichlet* bound =
        new chi_diffusion::BoundaryDirichlet();

      dsolver->boundaries[b] = bound;
    }//for bndry

    //================================= Redirect material lookup to use
    //                                  transport cross-sections
    dsolver->G  = 1;
    dsolver->gi = 0;

    //================================= Initialize solver, assemble matrix A
    //                                  but suppress solution
    bool verbose          = groupset->tgdsa_verbose;   //Disable normal info printing
    bool supress_assembly = false;   //Assemble the matrix
    bool supress_solver   = true;    //Suppress the solving
    dsolver->Initialize(verbose);
    dsolver->ExecuteS(supress_assembly,supress_solver);

    delta_phi_local.resize(0);
    delta_phi_local.shrink_to_fit();
  }//if wgdsa
}

//###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void CHI_NPTRANSPORT::AssembleTGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                                  double *ref_phi_old,
                                                  double *ref_phi_new)
{
  int gsi = groupset->groups[0]->id;
  int gss = groupset->groups.size();

  delta_phi_local.resize(local_dof_count,0.0);

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      chi_math::SparseMatrix& S = material_xs[xs_id]->transfer_matrix[0];

      for (int i=0; i<2; i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,0); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          double R_g = 0.0;
          int num_transfers = S.inds_rowI[gsi+g].size();
          for (int j=0; j<num_transfers; j++)
          {
            int gp = S.inds_rowI[gsi+g][j];

            if (not (gp >= (gsi+g+1)))
              continue;

            double delta_phi = phi_new_mapped[gp] - phi_old_mapped[gp];

            R_g += S.rowI_colJ[gsi+g][j]*delta_phi;
          }
          delta_phi_local[index] += R_g;
        }//for g

      }//for dof
    }//slab
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYGON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      chi_math::SparseMatrix& S = material_xs[xs_id]->transfer_matrix[0];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,0); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          double R_g = 0.0;
          int num_transfers = S.inds_rowI[gsi+g].size();
          for (int j=0; j<num_transfers; j++)
          {
            int gp = S.inds_rowI[gsi+g][j];

            if (not (gp >= (gsi+g+1)))
              continue;

            double delta_phi = phi_new_mapped[gp] - phi_old_mapped[gp];

            R_g += S.rowI_colJ[gsi+g][j]*delta_phi;
          }
          delta_phi_local[index] += R_g;
        }//for g

      }//for dof
    }//polygon
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYHEDRON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      chi_math::SparseMatrix& S = material_xs[xs_id]->transfer_matrix[0];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,0); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          double R_g = 0.0;
          int num_transfers = S.inds_rowI[gsi+g].size();
          for (int j=0; j<num_transfers; j++)
          {
            int gp = S.inds_rowI[gsi+g][j];

            if (not (gp >= (gsi+g+1)))
              continue;

            double delta_phi = phi_new_mapped[gp] - phi_old_mapped[gp];

            R_g += S.rowI_colJ[gsi+g][j]*delta_phi;
          }
          delta_phi_local[index] += R_g;
        }//for g

      }//for dof
    }//polyhedron

  }//for cell

}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void CHI_NPTRANSPORT::DisAssembleTGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                                     double *ref_phi_new)
{
  int gsi = groupset->groups[0]->id;
  int gss = groupset->groups.size();

  chi_diffusion::Solver* tgdsa_solver =
    (chi_diffusion::Solver*)groupset->tgdsa_solver;

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      std::vector<double>& xi_g = material_xs[xs_id]->xi_Jfull_g;


      for (int i=0; i<2; i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += tgdsa_solver->pwld_phi_local[index]*xi_g[gsi+g];

      }//for dof
    }//slab
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYGON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      std::vector<double>& xi_g = material_xs[xs_id]->xi_Jfull_g;


      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += tgdsa_solver->pwld_phi_local[index]*xi_g[gsi+g];

      }//for dof
    }//polygon
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYHEDRON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      std::vector<double>& xi_g = material_xs[xs_id]->xi_Jfull_g;


      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += tgdsa_solver->pwld_phi_local[index]*xi_g[gsi+g];

      }//for dof
    }//polyhedron


  }//for cell

  delta_phi_local.resize(0);
  delta_phi_local.shrink_to_fit();
}