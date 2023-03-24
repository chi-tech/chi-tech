#include "A_LBSSolver/lbs_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "A_LBSSolver/IterativeMethods/wgs_context.h"

#include "A_LBSSolver/Acceleration/diffusion_mip.h"



#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <iomanip>

namespace lbs
{

void PowerIterationKEigen1(LBSSolver& lbs_solver,
                          double tolerance,
                          int max_iterations,
                          double& k_eff)
{
  const std::string fname = "lbs::PowerIterationKEigen";

  for (auto& wgs_solver : lbs_solver.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context =
      std::dynamic_pointer_cast<lbs::WGSContext<Mat,Vec,KSP>>(context);

    if (not wgs_context) throw std::logic_error(fname + ": Cast failed.");

    wgs_context->lhs_src_scope_ = APPLY_WGS_SCATTER_SOURCES;
    wgs_context->rhs_src_scope_ = APPLY_AGS_SCATTER_SOURCES |
                                  APPLY_FIXED_SOURCES;
  }

  if (lbs_solver.Groupsets().size() > 1)
    throw std::invalid_argument(fname + ": The \"power1\" method can only "
                                        "be used with a single groupset.");

  auto& basic_options = lbs_solver.GetBasicOptions();
  auto& q_moments_local = lbs_solver.QMomentsLocal();
  auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& phi_new_local = lbs_solver.PhiNewLocal();
  auto primary_ags_solver = lbs_solver.GetPrimaryAGSSolver();
  auto& groupsets = lbs_solver.Groupsets();
  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();

  auto& front_gs = groupsets.front();
  auto& front_wgs_solver = lbs_solver.GetWGSSolvers()[front_gs.id_];
  auto frons_wgs_context = std::dynamic_pointer_cast<lbs::WGSContext<Mat,Vec,KSP>>(
    front_wgs_solver->GetContext());

  front_gs.apply_wgdsa_ = true;
  front_gs.wgdsa_tol_ = basic_options("PISA_MIP_L_ABS_TOL").FloatValue();
  front_gs.wgdsa_max_iters_ =
    static_cast<int>(basic_options("PISA_MIP_L_MAX_ITS").IntegerValue());
  lbs_solver.InitWGDSA(front_gs, true);
  front_gs.apply_wgdsa_ = false;

  bool pisa_verbose_level =
    static_cast<int>(basic_options("PISA_VERBOSE_LEVEL").IntegerValue());
  double pisa_pi_k_tol =
    basic_options("PISA_PI_K_TOL").FloatValue();
  int pisa_pi_max_its =
    static_cast<int>(basic_options("PISA_PI_MAX_ITS").IntegerValue());

  auto& diff_solver = *front_gs.wgdsa_solver_;

  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  /**Lambda for the creation of fission sources.*/
  auto SetLBSFissionSource = [&active_set_source_function,&front_gs]
    (const VecDbl& input, VecDbl& output)
  {
    chi_math::Set(output, 0.0);
    active_set_source_function(front_gs, output,
                               input,
                               APPLY_AGS_FISSION_SOURCES |
                               APPLY_WGS_FISSION_SOURCES);
  };

  //================================================== Start power iterations
  primary_ags_solver->SetVerbosity(lbs_solver.Options().verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    auto phi_l = lbs_solver.WGDSACopyOnlyPhi0(front_gs, phi_old_local);

    SetLBSFissionSource(/*input*/phi_old_local, /*output*/q_moments_local);
    chi_math::Scale(q_moments_local, 1.0/k_eff);

    auto Sffull = q_moments_local;
    auto Sf = lbs_solver.WGDSACopyOnlyPhi0(front_gs, q_moments_local);

    //====================================== This solves the inners for transport
    // produces phi at l+1/2
    primary_ags_solver->Setup();
    primary_ags_solver->Solve();

    //lph_i = l + 1/2,i
    auto phi_lph_i = lbs_solver.WGDSACopyOnlyPhi0(front_gs, phi_new_local);

    //lph_ip1 = l + 1/2, i+1
    q_moments_local = Sffull;
    active_set_source_function(front_gs, q_moments_local,
                               phi_new_local,
                               APPLY_AGS_SCATTER_SOURCES |
                               APPLY_WGS_SCATTER_SOURCES);

    frons_wgs_context->ApplyInverseTransportOperator(NO_FLAGS_SET);
    lbs_solver.GSScopedCopyPrimarySTLvectors(front_gs, phi_new_local, phi_old_local);

    auto phi_lph_ip1 = lbs_solver.WGDSACopyOnlyPhi0(front_gs, phi_new_local);

    //====================================== Power Iteration Acceleration solve
    auto epsilon_k = phi_l;
    chi_math::Set(epsilon_k, 0.0);
    auto epsilon_kp1 = epsilon_k;

    double mu_k = k_eff;
    double mu_kp1 = mu_k;

    for (size_t k=0; k<pisa_pi_max_its; ++k)
    {
      using namespace chi_math;

      auto eps_k_plus_phi_lph_ip1 = epsilon_k + phi_lph_ip1;
      lbs_solver.WGDSAProjectBackPhi0(front_gs, eps_k_plus_phi_lph_ip1, phi_old_local);

      SetLBSFissionSource(/*input*/phi_old_local, /*output*/q_moments_local);
      chi_math::Scale(q_moments_local, 1.0/mu_k);

      auto Sfaux = lbs_solver.WGDSACopyOnlyPhi0(front_gs, q_moments_local);

      Set(q_moments_local, 0.0);
      lbs_solver.WGDSAProjectBackPhi0(front_gs, phi_lph_ip1-phi_lph_i, phi_old_local);
      active_set_source_function(front_gs, q_moments_local,
                                 phi_old_local,
                                 APPLY_AGS_SCATTER_SOURCES |
                                 APPLY_WGS_SCATTER_SOURCES);

      auto Ss_res = lbs_solver.WGDSACopyOnlyPhi0(front_gs, q_moments_local);

      {
        Set(q_moments_local, 0.0);
        lbs_solver.WGDSAProjectBackPhi0(front_gs, epsilon_k, phi_old_local);
        active_set_source_function(front_gs, q_moments_local,
                                   phi_old_local,
                                   APPLY_AGS_SCATTER_SOURCES |
                                   APPLY_WGS_SCATTER_SOURCES |
                                   SUPPRESS_WG_SCATTER);

        auto Ss = lbs_solver.WGDSACopyOnlyPhi0(front_gs, q_moments_local);

        diff_solver.Assemble_b(Ss + Sfaux + Ss_res - Sf);
        diff_solver.Solve(epsilon_kp1);
      }

      lbs_solver.WGDSAProjectBackPhi0(front_gs, epsilon_k + phi_lph_ip1, phi_old_local);
      double pk = lbs_solver.ComputeFissionProduction(phi_old_local);

      lbs_solver.WGDSAProjectBackPhi0(front_gs, epsilon_kp1 + phi_lph_ip1, phi_old_local);
      double pkp1 = lbs_solver.ComputeFissionProduction(phi_old_local);

      mu_kp1 = pkp1/pk*mu_k;

      const double mu_change = std::fabs(1.0-mu_kp1/mu_k);
      if (pisa_verbose_level >= 1)
        chi::log.Log()
          << "PISA iteration " << k
          << " mu " << mu_kp1 << " mu change " << mu_change;

      if (mu_change < pisa_pi_k_tol) break;

      mu_k = mu_kp1;
      epsilon_k = epsilon_kp1;
    }//acceleration
    auto phi_lp1_temp = chi_math::operator+(epsilon_kp1, phi_lph_ip1);
    lbs_solver.WGDSAProjectBackPhi0(front_gs, phi_lp1_temp, phi_new_local);
    lbs_solver.GSScopedCopyPrimarySTLvectors(front_gs, phi_new_local, phi_old_local);

    const double production = lbs_solver.ComputeFissionProduction(phi_old_local);
    lbs_solver.ScalePhiVector(PhiSTLOption::PHI_OLD, mu_kp1/production);

    //======================================== Recompute k-eigenvalue
    k_eff = mu_kp1;
    double reactivity = (k_eff - 1.0) / k_eff;

    //======================================== Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    //======================================== Print iteration summary
    if (lbs_solver.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info
        << chi::program_timer.GetTimeString() << " "
        << "  Iteration " << std::setw(5) << nit
        << "  k_eff " << std::setw(11) << std::setprecision(7) << k_eff
        << "  k_eff change " << std::setw(12) << k_eff_change
        << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log()
    << "        Final k-eigenvalue    :        "
    << std::setprecision(7) << k_eff;
  chi::log.Log()
    << "        Final change          :        "
    << std::setprecision(6) << k_eff_change;
  chi::log.Log() << "\n";
}

}//namespace lbs