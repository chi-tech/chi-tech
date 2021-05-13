#include "lbs_linear_boltzmann_solver.h"

#include "../DiffusionSolver/Solver/diffusion_solver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initializes the Within-Group DSA solver. */
void LinearBoltzmann::Solver::InitTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa)
  {
    chi_math::UnknownManager scalar_uk_man;
    scalar_uk_man.AddUnknown(chi_math::UnknownType::SCALAR);

    //================================= Initialize field function
    delta_phi_local.resize(local_node_count, 0.0);
    int g = 0;
    int m = 0;
    std::string text_name = std::string("Sum_Sigma_s_DeltaPhi_g") +
                            std::to_string(g) +
                            std::string("_m") + std::to_string(m);

    auto deltaphi_ff = std::make_shared<chi_physics::FieldFunction>(
      text_name,                                    //Text name
      discretization,                               //Spatial Discretization
      &delta_phi_local,                             //Data vector
      scalar_uk_man);                               //Unknown manager

    chi_physics_handler.fieldfunc_stack.push_back(deltaphi_ff);
    field_functions.push_back(deltaphi_ff);

    //================================= Set diffusion solver
    std::string solver_name = std::string("TGDSA");
    solver_name += std::string("[g=")
                 + std::to_string(groupset.groups.front().id)
                 + std::string("-")
                 + std::to_string(groupset.groups.back().id)
                 + std::string("]");
    auto dsolver = new chi_diffusion::Solver(solver_name);
    groupset.tgdsa_solver = dsolver;

    dsolver->regions.push_back(this->regions.back());
    dsolver->discretization = discretization;
    dsolver->fem_method = PWLD_MIP;
    dsolver->residual_tolerance = groupset.tgdsa_tol;
    dsolver->max_iters          = groupset.tgdsa_max_iters;
    dsolver->options_string     = groupset.tgdsa_string;
    if (groupset.apply_wgdsa)
      dsolver->material_mode = DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL;
    else
      dsolver->material_mode = DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART;
    dsolver->q_field = deltaphi_ff;

    //================================= Initialize boundaries
    if (not dsolver->common_items_initialized)
      dsolver->InitializeCommonItems();

    typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
    dsolver->boundaries.clear();
    for (auto& lbs_bndry : sweep_boundaries)
    {
      if (lbs_bndry->Type() == SwpBndryType::REFLECTING)
        dsolver->boundaries.push_back(new chi_diffusion::BoundaryReflecting());
      else
        dsolver->boundaries.push_back(new chi_diffusion::BoundaryDirichlet());
    }

    //================================= Redirect material lookup to use
    //                                  transport cross-sections
    dsolver->G  = 1;
    dsolver->gi = 0;

    //================================= Initialize solver, assemble matrix A
    //                                  but suppress solution
    bool verbose          = groupset.tgdsa_verbose;   //Disable normal info printing
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
void LinearBoltzmann::Solver::CleanUpTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa)
    delete groupset.tgdsa_solver;
}

//###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void LinearBoltzmann::Solver::AssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                                          double *ref_phi_old,
                                                          double *ref_phi_new)
{
  int gsi = groupset.groups[0].id;
  int gss = groupset.groups.size();

  delta_phi_local.resize(local_node_count, 0.0);

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    int c = cell.local_id;

    auto& transport_view = cell_transport_views[c];

    int xs_id = matid_to_xs_map[cell.material_id];
    chi_math::SparseMatrix& S = material_xs[xs_id]->transfer_matrix[0];

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      index++;
      int m = 0;
      int mapping = transport_view.MapDOF(i,m,0); //phi_new & old location gsi

      double* phi_old_mapped = &ref_phi_old[mapping];
      double* phi_new_mapped = &ref_phi_new[mapping];

      for (int g=0; g<gss; g++)
      {
        double R_g = 0.0;
        int num_transfers = S.rowI_indices[gsi + g].size();
        for (int j=0; j<num_transfers; j++)
        {
          int gp = S.rowI_indices[gsi + g][j];

          if (gp < gsi + g + 1)
            continue;

          double delta_phi = phi_new_mapped[gp] - phi_old_mapped[gp];

          R_g += S.rowI_values[gsi + g][j] * delta_phi;
        }
        delta_phi_local[index] += R_g;
      }//for g

    }//for dof
  }//for cell

}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void LinearBoltzmann::Solver::DisAssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                                             double *ref_phi_new)
{
  int gsi = groupset.groups[0].id;
  int gss = groupset.groups.size();

  auto tgdsa_solver = (chi_diffusion::Solver*)groupset.tgdsa_solver;

  int index = -1;
  for (const auto& cell : grid->local_cells)
  {
    int c = cell.local_id;

    auto& transport_view = cell_transport_views[c];

    int xs_id = matid_to_xs_map[cell.material_id];
    std::vector<double>& xi_g = material_xs[xs_id]->xi_Jfull_g;

    for (int i=0; i < cell.vertex_ids.size(); i++)
    {
      index++;
      int m=0;
      int mapping = transport_view.MapDOF(i,m,gsi); //phi_new & old location gsi

      double* phi_new_mapped = &ref_phi_new[mapping];

      for (int g=0; g<gss; g++)
        phi_new_mapped[g] += tgdsa_solver->pwld_phi_local[index]*xi_g[gsi+g];

    }//for dof
  }//for cell

  delta_phi_local.resize(0);
  delta_phi_local.shrink_to_fit();
}