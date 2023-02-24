#include "A_LBSSolver/lbs_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "A_LBSSolver/Acceleration/diffusion_mip.h"

//###################################################################
/**Initializes the Within-Group DSA solver. */
void lbs::LBSSolver::InitWGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_wgdsa_)
  {
    //=========================================== Make UnknownManager
    const size_t gs_G = groupset.groups_.size();
    chi_math::UnknownManager uk_man;
    uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, gs_G);

    //=========================================== Make boundary conditions
    typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
    typedef lbs::acceleration::BoundaryCondition BC;
    typedef lbs::acceleration::BCType BCType;

    std::map<uint64_t, BC> bcs;
    for (auto& [bid,lbs_bndry] : sweep_boundaries_)
    {
      if (lbs_bndry->Type() == SwpBndryType::REFLECTING)
        bcs[bid] = {BCType::ROBIN,{0.0,1.0,0.0}};
      else//dirichlet
        bcs[bid] = {BCType::DIRICHLET,{0.0,0.0,0.0}};
    }

    //=========================================== Make xs map
    typedef lbs::acceleration::Multigroup_D_and_sigR MGXS;
    typedef std::map<int, lbs::acceleration::Multigroup_D_and_sigR> MatID2XSMap;
    MatID2XSMap matid_2_mgxs_map;
    for (const auto& matid_xs_pair : matid_to_xs_map_)
    {
      const auto& mat_id = matid_xs_pair.first;
      const auto& xs     = matid_xs_pair.second;

      std::vector<double> Dg  (gs_G, 0.0);
      std::vector<double> sigR(gs_G, 0.0);

      size_t g = 0;
      for (size_t gprime = groupset.groups_.front().id_;
           gprime <= groupset.groups_.back().id_; ++gprime)
      {
        Dg[g]   = xs->diffusion_coeff[gprime];
        sigR[g] = xs->sigma_removal[gprime];
        ++g;
      }//for g

      matid_2_mgxs_map.insert(std::make_pair(mat_id, MGXS{Dg, sigR}));
    }

    //=========================================== Create solver
    const auto& sdm = *discretization_;

    auto solver =
      std::make_shared<acceleration::DiffusionMIPSolver>(
        std::string(TextName()+"_WGDSA"),
        *grid_ptr_, sdm,
        uk_man,
        bcs,
        matid_2_mgxs_map,
        unit_cell_matrices_,
        true); //verbosity

    solver->options.residual_tolerance        = groupset.wgdsa_tol_;
    solver->options.max_iters                 = groupset.wgdsa_max_iters_;
    solver->options.verbose                   = groupset.wgdsa_verbose_;
    solver->options.additional_options_string = groupset.wgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man),0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.wgdsa_solver_ = solver;
  }
}

//###################################################################
/**Cleans up memory consuming items. */
void lbs::LBSSolver::CleanUpWGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_wgdsa_) groupset.wgdsa_solver_ = nullptr;
}


//###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void lbs::LBSSolver::
  AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                              const std::vector<double>& phi_in,
                              std::vector<double>& delta_phi_local)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man  = flux_moments_uk_man_;

  const int    gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  delta_phi_local.clear();
  delta_phi_local.assign(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& sigma_s = matid_to_xs_map_[cell.material_id]->sigma_s_gtog;

    for (size_t i=0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t  phi_map = sdm.MapDOFLocal(cell, i,  phi_uk_man, 0, gsi);

      double* delta_phi_mapped    = &delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g=0; g<gss; g++)
      {
        delta_phi_mapped[g] = sigma_s[gsi+g] * phi_in_mapped[g];
      }//for g
    }//for node
  }//for cell
}

//###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void lbs::LBSSolver::
  DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                 const std::vector<double>& delta_phi_local,
                                 std::vector<double>& ref_phi_new)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man  = flux_moments_uk_man_;

  const int    gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i=0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t  phi_map = sdm.MapDOFLocal(cell, i,  phi_uk_man, 0, gsi);

      const double* delta_phi_mapped = &delta_phi_local[dphi_map];
            double* phi_new_mapped   = &ref_phi_new[phi_map];

      for (int g=0; g<gss; g++)
        phi_new_mapped[g] += delta_phi_mapped[g];
    }//for dof
  }//for cell
}