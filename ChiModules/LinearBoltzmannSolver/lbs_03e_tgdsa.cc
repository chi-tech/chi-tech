#include "lbs_linear_boltzmann_solver.h"

#include "Groupset/lbs_groupset.h"

#include "Acceleration/diffusion_mip.h"

//###################################################################
/**Initializes the Within-Group DSA solver. */
void lbs::SteadySolver::InitTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa)
  {
    //=========================================== Make UnknownManager
    const auto& uk_man = discretization->UNITARY_UNKNOWN_MANAGER;

    //=========================================== Make boundary conditions
    typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
    typedef lbs::acceleration::BoundaryCondition BC;
    typedef lbs::acceleration::BCType BCType;

    std::vector<BC> bcs;
    size_t bid = 0;
    for (auto& lbs_bndry : sweep_boundaries)
    {
      if (lbs_bndry->Type() == SwpBndryType::REFLECTING)
        bcs.push_back({BCType::ROBIN,{0.0,1.0,0.0}});
      else//dirichlet
        bcs.push_back({BCType::DIRICHLET,{0.0,0.0,0.0}});
      ++bid;
    }

    //=========================================== Make xs map
    typedef lbs::acceleration::Multigroup_D_and_sigR MGXs;
    typedef std::map<int, lbs::acceleration::Multigroup_D_and_sigR> MapMatID2XS;
    MapMatID2XS map_mat_id_2_mgxs;
    for (const auto& mat_id_xs_pair : matid_to_xs_map)
    {
      const auto& mat_id = mat_id_xs_pair.first;
      const auto& xs     = mat_id_xs_pair.second;

      std::vector<double> Dg  (1, 0.0);
      std::vector<double> sigR(1, 0.0);

      Dg[0]   = xs->D_jfull;
      sigR[0] = xs->sigma_a_jfull;

      map_mat_id_2_mgxs.insert(std::make_pair(mat_id,MGXs{Dg,sigR}));
    }

    //=========================================== Create solver
    const auto& sdm = *discretization;

    auto solver =
      std::make_shared<acceleration::DiffusionMIPSolver>(
        std::string(TextName()+"_TGDSA"),
        *grid,sdm,
        uk_man,
        bcs,
        map_mat_id_2_mgxs,
        unit_cell_matrices,
        true); //verbosity

    solver->options.residual_tolerance        = groupset.tgdsa_tol;
    solver->options.max_iters                 = groupset.tgdsa_max_iters;
//    solver->options.verbose                   = groupset.tgdsa_verbose;
    solver->options.additional_options_string = groupset.tgdsa_string;

    solver->Initialize();

    delta_phi_local.assign(sdm.GetNumLocalDOFs(uk_man),0.0);

    solver->AssembleAand_b2(delta_phi_local);

    delta_phi_local.resize(0);
    delta_phi_local.shrink_to_fit();

    groupset.tgdsa_solver = solver;
  }
}

//###################################################################
/**Cleans up memory consuming items. */
void lbs::SteadySolver::CleanUpTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa) groupset.tgdsa_solver = nullptr;
}

//###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void lbs::SteadySolver::AssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                                    const double *ref_phi_old,
                                                    const double *ref_phi_new)
{
  const auto& sdm = *discretization;
  const auto& dphi_uk_man = groupset.tgdsa_solver->UnknownStructure();
  const auto& phi_uk_man  = flux_moments_uk_man;

  const int    gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  delta_phi_local.assign(local_node_count, 0.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& S = matid_to_xs_map[cell.material_id]->transfer_matrices[0];

    for (size_t i=0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t  phi_map = sdm.MapDOFLocal(cell, i,  phi_uk_man, 0, gsi);

      const double* phi_old_mapped   = &ref_phi_old[phi_map];
      const double* phi_new_mapped   = &ref_phi_new[phi_map];

      for (size_t g=0; g<gss; g++)
      {
        delta_phi_local[dphi_map] = 0.0;
        double R_g = 0.0;
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(g))
          R_g += sigma_sm * (phi_new_mapped[gprime] - phi_old_mapped[gprime]);

        delta_phi_local[dphi_map] += R_g;
      }//for g
    }//for node
  }//for cell
}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void lbs::SteadySolver::DisAssembleTGDSADeltaPhiVector(LBSGroupset& groupset,
                                                       double *ref_phi_new)
{
  const auto& sdm = *discretization;
  const auto& dphi_uk_man = groupset.tgdsa_solver->UnknownStructure();
  const auto& phi_uk_man  = flux_moments_uk_man;

  const int    gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    auto& transport_view = cell_transport_views[cell.local_id];

    const auto& xi_g = matid_to_xs_map[cell.material_id]->xi_Jfull;

    for (size_t i=0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t  phi_map = sdm.MapDOFLocal(cell, i,  phi_uk_man, 0, gsi);

      const double delta_phi_mapped = delta_phi_local[dphi_map];
      double* phi_new_mapped   = &ref_phi_new[phi_map];

      for (int g=0; g<gss; g++)
        phi_new_mapped[g] += delta_phi_mapped*xi_g[gsi+g];
    }//for dof
  }//for cell

  delta_phi_local.resize(0);
  delta_phi_local.shrink_to_fit();
}

//###################################################################
/**Executes Two-Grid Acceleration. This involves assembling the system RHS,
 * solving the system and finally adding the solution to the scalar flux.*/
void lbs::SteadySolver::
  ExecuteTGDSA(LBSGroupset &groupset,
               const std::vector<double>& ref_phi_old,
               std::vector<double>& ref_phi_new)
{
  AssembleTGDSADeltaPhiVector(groupset,
                              ref_phi_old.data(),
                              ref_phi_new.data());
  groupset.tgdsa_solver->Assemble_b2(delta_phi_local);
  groupset.tgdsa_solver->Solve(delta_phi_local);
  DisAssembleTGDSADeltaPhiVector(groupset, ref_phi_new.data());
}