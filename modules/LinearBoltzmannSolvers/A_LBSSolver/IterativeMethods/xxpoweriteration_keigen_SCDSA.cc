#include "xxpoweriteration_keigen_SCDSA.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "A_LBSSolver/Acceleration/diffusion_mip.h"
#include "A_LBSSolver/Acceleration/diffusion_PWLC.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include <iomanip>

namespace lbs
{

RegisterChiObject(lbs, XXPowerIterationKEigenSCDSA);

chi_objects::InputParameters XXPowerIterationKEigenSCDSA::GetInputParameters()
{
  chi_objects::InputParameters params =
    XXPowerIterationKEigen::GetInputParameters();

  params.AddOptionalParameter(
    "accel_pi_max_its",
    50,
    "Maximum allowable iterations for the acceleration scheme's inner "
    "power iterations");

  params.AddOptionalParameter(
    "accel_pi_k_tol",
    1.0e-10,
    "K-eigenvalue tolerance for the acceleration scheme's inner "
    "power iterations");

  params.AddOptionalParameter(
    "accel_pi_verbose",
    false,
    "Flag, if set will result in verbose output from the acceleration "
    "scheme");

  params.AddOptionalParameter(
    "diff_accel_diffusion_l_abs_tol",
    1.0e-10,
    "Absolute residual tolerance to use for the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_max_iters",
    100,
    "Maximum allowable iterations for the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_verbose",
    false,
    "Flag, if set will enable verbose output of the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_petsc_options",
    std::string("ssss"),
    "Additional PETSc options for the diffusion accelerator");

  return params;
}

XXPowerIterationKEigenSCDSA::XXPowerIterationKEigenSCDSA(
  const chi_objects::InputParameters& params)
  : XXPowerIterationKEigen(params),
    accel_pi_max_its_(params.GetParamValue<int>("accel_pi_max_its")),
    accel_pi_k_tol_(params.GetParamValue<double>("accel_pi_k_tol")),
    accel_pi_verbose_(params.GetParamValue<bool>("accel_pi_verbose"))
{
  //=========================================== Make UnknownManager
  const size_t num_gs_groups = front_gs_.groups_.size();
  chi_math::UnknownManager uk_man;
  uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_gs_groups);

  //=========================================== Make boundary conditions
  auto bcs = acceleration::TranslateBCs(lbs_solver_.SweepBoundaries(),
                                        /*vaccum_bcs_are_dirichlet=*/false);

  //=========================================== Make xs map
  auto matid_2_mgxs_map =
    acceleration::PackGroupsetXS(lbs_solver_.GetMatID2XSMap(),
                                 front_gs_.groups_.front().id_,
                                 front_gs_.groups_.back().id_);

  //=========================================== Create solver
  const auto& sdm = lbs_solver_.SpatialDiscretization();
  const auto& unit_cell_matrices = lbs_solver_.GetUnitCellMatrices();

  diffusion_solver_ = std::make_shared<acceleration::DiffusionMIPSolver>(
    std::string(TextName() + "_WGDSA"),
    sdm,
    uk_man,
    bcs,
    matid_2_mgxs_map,
    unit_cell_matrices,
    true); // verbosity
  chi_objects::ParameterBlock block;

  {
    typedef const std::string cstr;
    cstr l_abs_tol = "diff_accel_diffusion_l_abs_tol";
    cstr max_iters = "diff_accel_diffusion_max_iters";
    cstr verbose = "diff_accel_diffusion_verbose";
    cstr petsc_options = "diff_accel_diffusion_petsc_options";

    auto& ds = diffusion_solver_;

    ds->options.residual_tolerance = params.GetParamValue<double>(l_abs_tol);
    ds->options.max_iters = params.GetParamValue<int>(max_iters);
    ds->options.verbose = params.GetParamValue<bool>(verbose);
    if (not params.GetParamValue<std::string>(petsc_options).empty())
      ds->options.additional_options_string =
        params.GetParamValue<std::string>(petsc_options);
  }

  diffusion_solver_->Initialize();

  std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);
  diffusion_solver_->AssembleAand_b(dummy_rhs);
}

// ##################################################################
/**Executes the scheme.*/
void XXPowerIterationKEigenSCDSA::Execute()
{
  auto phi_temp = phi_old_local_;

  /**Lambda for the creation of scattering sources but the
   * input vector is only the zeroth moment*/
  auto SetLBSScatterSourcePhi0 =
    [this, &phi_temp](
      const VecDbl& input, const bool additive, const bool suppress_wgs = false)
  {
    ProjectBackPhi0(front_gs_, input, phi_temp);
    SetLBSScatterSource(/*in*/ phi_temp, additive, suppress_wgs);
  };

  using namespace chi_math;

  k_eff_ = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iters_)
  {
    auto phi0_l = CopyOnlyPhi0(front_gs_, phi_old_local_);

    //================================= Set the fission source
    SetLBSFissionSource(phi_old_local_, /*additive=*/false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    auto Sf_all_moments = q_moments_local_;
    auto Sf = CopyOnlyPhi0(front_gs_, q_moments_local_);

    //================================= This solves the inners for transport
    primary_ags_solver_.Setup();
    primary_ags_solver_.Solve();

    // lph_i = l + 1/2,i
    auto phi0_lph_i = CopyOnlyPhi0(front_gs_, phi_new_local_);

    // Now we produce lph_ip1 = l + 1/2, i+1
    q_moments_local_ = Sf_all_moments; // Restore 1/k F phi_l
    SetLBSScatterSource(phi_new_local_, /*additive=*/true);

    front_wgs_context_->ApplyInverseTransportOperator(NO_FLAGS_SET); // Sweep

    auto phi0_lph_ip1 = CopyOnlyPhi0(front_gs_, phi_new_local_);

    //====================================== Power Iteration Acceleration solve
    VecDbl epsilon_k(phi0_l.size(), 0.0);
    auto epsilon_kp1 = epsilon_k;

    double lambda_k = k_eff_;
    double lambda_kp1 = lambda_k;

    for (size_t k = 0; k < accel_pi_max_its_; ++k)
    {
      ProjectBackPhi0(front_gs_,
                      /*in*/ epsilon_k + phi0_lph_ip1,
                      /*out*/ phi_temp);

      double production_k = lbs_solver_.ComputeFissionProduction(phi_temp);

      SetLBSFissionSource(phi_temp, /*additive=*/false);
      Scale(q_moments_local_, 1.0 / lambda_k);

      auto Sfaux = CopyOnlyPhi0(front_gs_, q_moments_local_);

      SetLBSScatterSourcePhi0(phi0_lph_ip1 - phi0_lph_i, /*additive=*/false);

      auto Ss_res = CopyOnlyPhi0(front_gs_, q_moments_local_);

      // Inner iterations seems extremely wasteful therefore I
      // am leaving this at 1 iteration here for further investigation.
      for (int i = 0; i < 1; ++i)
      {
        SetLBSScatterSourcePhi0(epsilon_k,
                                /*additive=*/false,
                                /*suppress_wgs=*/true);

        auto Ss = CopyOnlyPhi0(front_gs_, q_moments_local_);

        // Solve the diffusion system
        diffusion_solver_->Assemble_b(Ss + Sfaux + Ss_res - Sf);
        diffusion_solver_->Solve(epsilon_kp1, /*use_initial_guess=*/true);

        epsilon_k = epsilon_kp1;
      }

      ProjectBackPhi0(front_gs_,
                      /*in*/ epsilon_kp1 + phi0_lph_ip1,
                      /*out*/ phi_old_local_);
      double production_kp1 =
        lbs_solver_.ComputeFissionProduction(phi_old_local_);

      lambda_kp1 = production_kp1 / (production_k / lambda_k);

      const double lambda_change = std::fabs(1.0 - lambda_kp1 / lambda_k);
      if (accel_pi_verbose_ >= 1)
        chi::log.Log() << "PISA iteration " << k << " lambda " << lambda_kp1
                       << " lambda change " << lambda_change;

      if (lambda_change < accel_pi_k_tol_) break;

      lambda_k = lambda_kp1;
      epsilon_k = epsilon_kp1;
    } // acceleration

    ProjectBackPhi0(front_gs_,
                    /*in*/ epsilon_kp1 + phi0_lph_ip1,
                    /*out*/ phi_new_local_);
    lbs_solver_.GSScopedCopyPrimarySTLvectors(front_gs_,
                                              /*in*/ phi_new_local_,
                                              /*out*/ phi_old_local_);

    const double production =
      lbs_solver_.ComputeFissionProduction(phi_old_local_);
    lbs_solver_.ScalePhiVector(PhiSTLOption::PHI_OLD, lambda_kp1 / production);

    //================================= Recompute k-eigenvalue
    k_eff_ = lambda_kp1;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    //================================= Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    nit += 1;

    if (k_eff_change < std::max(k_tolerance_, 1.0e-12)) converged = true;

    //================================= Print iteration summary
    if (lbs_solver_.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << chi::program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff "
                  << std::setw(11) << std::setprecision(7) << k_eff_
                  << "  k_eff change " << std::setw(12) << k_eff_change
                  << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  } // for k iterations

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log() << "        Final k-eigenvalue    :        "
                 << std::setprecision(7) << k_eff_;
  chi::log.Log() << "        Final change          :        "
                 << std::setprecision(6) << k_eff_change << " (num_TrOps:"
                 << front_wgs_context_->counter_applications_of_inv_op_ << ")"
                 << "\n";
  chi::log.Log() << "\n";

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    chi_math::Scale(lbs_solver_.PrecursorsNewLocal(), 1.0 / k_eff_);
  }

  lbs_solver_.UpdateFieldFunctions();

  chi::log.Log()
    << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

// ##################################################################
/**Copies only the scalar moments from an lbs primary flux moments
* vector.*/
std::vector<double>
XXPowerIterationKEigenSCDSA::CopyOnlyPhi0(const LBSGroupset& groupset,
                                          const std::vector<double>& phi_in)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  VecDbl output_phi_local(diff_sdm.GetNumLocalDOFs(diff_uk_man), 0.0);

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[diff_phi_map];
      const double* phi_in_mapped = &phi_in[lbs_phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

// ##################################################################
/**Copies back only the scalar moments to a lbs primary flux vector.*/
void XXPowerIterationKEigenSCDSA::ProjectBackPhi0(
  const LBSGroupset& groupset,
  const std::vector<double>& input,
  std::vector<double>& output)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[diff_phi_map];
      double* output_mapped = &output[lbs_phi_map];

      for (int g = 0; g < gss; g++)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

} // namespace lbs