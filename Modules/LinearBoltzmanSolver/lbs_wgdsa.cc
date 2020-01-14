#include "lbs_linear_boltzman_solver.h"

#include <ChiMesh/Cell/cell.h>

#include "../DiffusionSolver/Solver/diffusion_solver.h"
#include "../DiffusionSolver/Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../DiffusionSolver/Boundaries/chi_diffusion_bndry_reflecting.h"
#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Initializes the Within-Group DSA solver. */
void LinearBoltzman::Solver::InitWGDSA(LBSGroupset *groupset)
{
  if (groupset->apply_wgdsa)
  {
    //================================= Initialize field function
    delta_phi_local.resize(local_dof_count*groupset->groups.size(),0.0);
    int g = 0;
    int m = 0;
    auto deltaphi_ff = new chi_physics::FieldFunction;
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
    auto dsolver = new chi_diffusion::Solver(solver_name);
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

    typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
    dsolver->boundaries.clear();
    for (auto lbs_bndry : sweep_boundaries)
    {
      if (lbs_bndry->Type() == SwpBndryType::REFLECTING)
        dsolver->boundaries.push_back(new chi_diffusion::BoundaryReflecting());
      else
        dsolver->boundaries.push_back(new chi_diffusion::BoundaryDirichlet());
    }

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
/**Cleans up memory consuming items. */
void LinearBoltzman::Solver::CleanUpWGDSA(LBSGroupset *groupset)
{
  if (groupset->apply_wgdsa)
    delete groupset->wgdsa_solver;
}

//###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void LinearBoltzman::Solver::AssembleWGDSADeltaPhiVector(LBSGroupset *groupset,
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

    auto transport_view =
      (LinearBoltzman::CellViewFull*)cell_transport_views[c];

    int xs_id = matid_to_xs_map[cell->material_id];
    std::vector<double>& sigma_s = material_xs[xs_id]->sigma_s_gtog;

    for (int i=0; i < cell->vertex_ids.size(); i++)
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
  }//for cell

}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void LinearBoltzman::Solver::DisAssembleWGDSADeltaPhiVector(LBSGroupset *groupset,
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

    auto transport_view =
      (LinearBoltzman::CellViewFull*)cell_transport_views[c];

    for (int i=0; i < cell->vertex_ids.size(); i++)
    {
      index++;
      int m=0;
      int mapping = transport_view->MapDOF(i,m,gsi); //phi_new & old location gsi

      double* phi_new_mapped = &ref_phi_new[mapping];

      for (int g=0; g<gss; g++)
        phi_new_mapped[g] += wgdsa_solver->pwld_phi_local[index*gss+g];

    }//for dof
  }//for cell

  delta_phi_local.resize(0);
  delta_phi_local.shrink_to_fit();
}