#include "lbs_linear_boltzman_solver.h"

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
void LinearBoltzmanSolver::InitWGDSA(NPT_GROUPSET *groupset)
{
  if (groupset->apply_wgdsa)
  {
    //================================= Initialize field function
    delta_phi_local.resize(local_dof_count*groupset->groups.size(),0.0);
    int g = 0;
    int m = 0;
    chi_physics::FieldFunction* deltaphi_ff =
      new chi_physics::FieldFunction;
    deltaphi_ff->text_name = std::string("Sigma_s_DeltaPhi_g") +
                          std::to_string(g) +
                          std::string("_m") + std::to_string(m);
    deltaphi_ff->grid = grid;
    deltaphi_ff->spatial_discretization = discretization;
    deltaphi_ff->id = chi_physics_handler.fieldfunc_stack.size();

    deltaphi_ff->type = FF_SDM_PWLD;
    deltaphi_ff->num_grps = groupset->groups.size();
    deltaphi_ff->num_moms = 1;
    deltaphi_ff->grp = g;
    deltaphi_ff->mom = m;
    deltaphi_ff->field_vector_local = &delta_phi_local;

    groupset->wgdsa_cell_dof_array_address.resize(
      local_cell_dof_array_address.size(),0);
    for (int c=0; c<local_cell_dof_array_address.size(); c++)
    {
      groupset->wgdsa_cell_dof_array_address[c] =
        local_cell_dof_array_address[c]*groupset->groups.size();
    }

    deltaphi_ff->local_cell_dof_array_address =
      &groupset->wgdsa_cell_dof_array_address;

    chi_physics_handler.fieldfunc_stack.push_back(deltaphi_ff);
    field_functions.push_back(deltaphi_ff);

    //================================= Set diffusion solver
    std::string solver_name = std::string("WGDSA");
    chi_diffusion::Solver* dsolver = new chi_diffusion::Solver(solver_name);
    groupset->wgdsa_solver = dsolver;

    dsolver->regions.push_back(this->regions.back());
    dsolver->discretization = discretization;
    dsolver->fem_method = PWLD_MIP_GAGG;
    dsolver->residual_tolerance = groupset->wgdsa_tol;
    dsolver->max_iters          = groupset->wgdsa_max_iters;
    dsolver->options_string     = groupset->wgdsa_string;
    dsolver->material_mode = DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF;
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
    dsolver->G  = groupset->groups.size();
    dsolver->gi = groupset->groups.front()->id;

    //================================= Initialize solver, assemble matrix A
    //                                  but suppress solution
    bool verbose          = groupset->wgdsa_verbose;   //Disable normal info printing
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
void LinearBoltzmanSolver::AssembleWGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                                  double *ref_phi_old,
                                                  double *ref_phi_new)
{
  int gsi = groupset->groups[0]->id;
  int gss = groupset->groups.size();

  delta_phi_local.resize(local_dof_count*gss,0.0);

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      int xs_id = matid_to_xs_map[cell->material_id];
      std::vector<double>& sigma_s = material_xs[xs_id]->sigma_s_gtog;

      for (int i=0; i<2; i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          delta_phi_local[index*gss+g] =
            phi_new_mapped[g] - phi_old_mapped[g];
          delta_phi_local[index*gss+g] *= sigma_s[gsi+g];
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
      std::vector<double>& sigma_s = material_xs[xs_id]->sigma_s_gtog;

      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          delta_phi_local[index*gss+g] =
            phi_new_mapped[g] - phi_old_mapped[g];
          delta_phi_local[index*gss+g] *= sigma_s[gsi+g];
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
      std::vector<double>& sigma_s = material_xs[xs_id]->sigma_s_gtog;

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        index++;
        int m = 0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_old_mapped = &ref_phi_old[mapping];
        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
        {
          delta_phi_local[index*gss+g] =
            phi_new_mapped[g] - phi_old_mapped[g];
          delta_phi_local[index*gss+g] *= sigma_s[gsi+g];
        }//for g
      }//for dof
    }//polyhedron
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encounted in call to "
           "AssembleWGDSADeltaPhiVector.";
      exit(EXIT_FAILURE);
    }

  }//for cell

}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void LinearBoltzmanSolver::DisAssembleWGDSADeltaPhiVector(NPT_GROUPSET *groupset,
                                                     double *ref_phi_new)
{
  int gsi = groupset->groups[0]->id;
  int gss = groupset->groups.size();

  chi_diffusion::Solver* wgdsa_solver =
    (chi_diffusion::Solver*)groupset->wgdsa_solver;

  int index = -1;
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell        = grid->cells[cell_g_index];

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<2; i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += wgdsa_solver->pwld_phi_local[index*gss+g];

      }//for dof
    }//slab
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYGON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += wgdsa_solver->pwld_phi_local[index*gss+g];

      }//for dof
    }//polygon
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& POLYHEDRON
    else if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      NPT_CELLVIEW_FULL* transport_view =
        (NPT_CELLVIEW_FULL*)cell_transport_views[c];

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        index++;
        int m=0;
        int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

        double* phi_new_mapped = &ref_phi_new[mapping];

        for (int g=0; g<gss; g++)
          phi_new_mapped[g] += wgdsa_solver->pwld_phi_local[index*gss+g];

      }//for dof
    }//polyhedron
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encounted in call to "
           "DisAssembleWGDSADeltaPhiVector.";
      exit(EXIT_FAILURE);
    }

  }//for cell

  delta_phi_local.resize(0);
  delta_phi_local.shrink_to_fit();
}