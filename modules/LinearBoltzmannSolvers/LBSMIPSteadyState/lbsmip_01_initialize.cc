#include "lbsmip_steady_solver.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"
#include "IterativeOperations/mip_wgs_context.h"
#include "LBSSteadyState/IterativeMethods/wgs_linear_solver.h"

/**Initializing.*/
void lbs::MIPSteadyStateSolver::Initialize()
{
  options_.scattering_order = 0; //overwrite any setting otherwise
  SteadyStateSolver::Initialize();
}

/**Initializes Within-GroupSet solvers.*/
void lbs::MIPSteadyStateSolver::InitializeWGSSolvers()
{
  //============================================= Initialize groupset solvers
  gs_mip_solvers_.assign(groupsets_.size(), nullptr);
  const size_t num_groupsets = groupsets_.size();
  for (size_t gs=0; gs<num_groupsets; ++gs)
  {
    const auto& groupset = groupsets_[gs];

    //=========================================== Make UnknownManager
    const size_t gs_G = groupset.groups.size();
    chi_math::UnknownManager uk_man;
    uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, gs_G);

    //=========================================== Make boundary conditions
    typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
    typedef lbs::acceleration::BoundaryCondition BC;
    typedef lbs::acceleration::BCType BCType;

    std::map<uint64_t, BC> bcs;
    for (auto& [bid, lbs_bndry] : sweep_boundaries_)
    {
      if (lbs_bndry->Type() == SwpBndryType::REFLECTING)
        bcs[bid] = {BCType::ROBIN,{0.0,1.0,0.0}};
      else if (lbs_bndry->Type() == SwpBndryType::VACUUM)
        bcs[bid] = {BCType::ROBIN,{0.25,0.5}};
    }//for sweep-boundary

    //=========================================== Make xs map
    typedef lbs::acceleration::Multigroup_D_and_sigR MGXs;
    typedef std::map<int, lbs::acceleration::Multigroup_D_and_sigR> MapMatID2XS;
    MapMatID2XS map_mat_id_2_mgxs;
    for (const auto& mat_id_xs_pair : matid_to_xs_map_)
    {
      const auto& mat_id = mat_id_xs_pair.first;
      const auto& xs     = mat_id_xs_pair.second;

      std::vector<double> Dg  (gs_G, 0.0);
      std::vector<double> sigR(gs_G, 0.0);

      size_t g = 0;
      for (size_t gprime=groupset.groups.front().id;
           gprime<=groupset.groups.back().id; ++gprime)
      {
        Dg[g]   = xs->diffusion_coeff[gprime];
        sigR[g] = xs->sigma_removal[gprime];
        ++g;
      }//for g

      map_mat_id_2_mgxs.insert(std::make_pair(mat_id,MGXs{Dg,sigR}));
    }

    //=========================================== Create solver
    const auto& sdm = *discretization_;

    auto solver =
      std::make_shared<acceleration::DiffusionMIPSolver>(
        std::string(TextName()+"_WGSolver"),
        *grid_ptr_, sdm,
        uk_man,
        bcs,
        map_mat_id_2_mgxs,
        unit_cell_matrices_,
        true); //verbosity

    solver->options.residual_tolerance        = groupset.wgdsa_tol;
    solver->options.max_iters                 = groupset.wgdsa_max_iters;
    solver->options.verbose                   = groupset.wgdsa_verbose;
    solver->options.additional_options_string = groupset.wgdsa_string;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man),0.0);

    solver->AssembleAand_b(dummy_rhs);

    gs_mip_solvers_[gs] = solver;
  }//for groupset

  wgs_solvers_.clear(); //this is required
  for (auto& groupset : groupsets_)
  {

    auto mip_wgs_context_ptr =
    std::make_shared<MIPWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES |
        SUPPRESS_WG_SCATTER,                                    //lhs_scope
        APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
        APPLY_AGS_FISSION_SOURCES,                             //rhs_scope
        options_.verbose_inner_iterations);

    auto wgs_solver =
      std::make_shared<WGSLinearSolver<Mat,Vec,KSP>>(mip_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  }//for groupset

}