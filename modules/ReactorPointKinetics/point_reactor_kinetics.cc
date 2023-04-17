#include "point_reactor_kinetics.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <numeric>

namespace prk
{

RegisterChiObject(prk, TransientSolver);

/**Sets input parameters.*/
chi_objects::InputParameters TransientSolver::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "rpk_TransientSolver");

  std::vector<double> default_lambdas = {
    0.0124, 0.0304, 0.111, 0.301, 1.14, 3.01};
  std::vector<double> default_betas = {
    0.00021, 0.00142, 0.00127, 0.00257, 0.00075, 0.00027};

  params.AddOptionalParameterBlock(
    "precursor_lambdas", default_lambdas, "An array of decay constants");
  params.AddOptionalParameterBlock(
    "precursor_betas",
    default_betas,
    "An array of fractional delayed neutron fractions");
  params.AddOptionalParameter(
    "gen_time", 1.0e-5, "Neutron generation time [s]");
  params.AddOptionalParameter("initial_rho", 0.0, "Initial reactivity [$]");
  params.AddOptionalParameter(
    "initial_source", 1.0, "Initial source strength [/s]");
  params.AddOptionalParameter(
    "initial_population", 1.0, "Initial neutron population");
  params.AddOptionalParameter("dt", 0.01, "Default timestep size [s]");

  params.AddOptionalParameter(
    "time_integration", "implicit_euler", "Time integration scheme to use");

  using namespace chi_data_types;
  auto time_intgl_list = AllowableRangeList::New(
    {"explicit_euler", "implicit_euler", "crank_nicolson"});

  params.ConstrainParameterRange("time_integration",
                                 std::move(time_intgl_list));

  params.ConstrainParameterRange("gen_time",
                                 AllowableRangeLowLimit::New(1.0e-12));
  params.ConstrainParameterRange(
    "dt", AllowableRangeLowHighLimit::New(1.0e-12, 100.0));
  params.ConstrainParameterRange("initial_source",
                                 AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("initial_population",
                                 AllowableRangeLowLimit::New(0.0));
  return params;
}

/**Constructor.*/
TransientSolver::TransientSolver(const chi_objects::InputParameters& params)
  : chi_physics::Solver(params.GetParamValue<std::string>("name")),
    lambdas_(params.GetParamVectorValue<double>("precursor_lambdas")),
    betas_(params.GetParamVectorValue<double>("precursor_betas")),
    gen_time_(params.GetParamValue<double>("gen_time")),
    rho_(params.GetParamValue<double>("initial_rho")),
    source_strength_(params.GetParamValue<double>("initial_source")),
    num_precursors_(lambdas_.size()),
    dt_(params.GetParamValue<double>("dt"))
{
  chi::log.Log() << "Created solver " << TextName();
  {
    std::stringstream outstr;
    outstr << "lambdas = ";
    for (double val : lambdas_)
      outstr << val << " ";
    chi::log.Log() << outstr.str();
  }
  {
    std::stringstream outstr;
    outstr << "betas = ";
    for (double val : betas_)
      outstr << val << " ";
    chi::log.Log() << outstr.str();
  }
}

/**Initialize function.*/
void TransientSolver::Initialize()
{
  // Check size
  if (lambdas_.size() != betas_.size())
    throw std::logic_error(TextName() +
                           ": Number of precursors cannot be "
                           "deduced from precursor data because "
                           "the data lists are of different size.");

  beta_ = std::accumulate(betas_.begin(), betas_.end(), /*init_val=*/0.0);

  // ============================= Initializing linalg items
  const auto& J = num_precursors_;
  A_ = chi_math::DynamicMatrix<double>(J + 1, J + 1, 0.0);
  I_ = A_;
  I_.SetDiagonal(1.0);

  x_t_ = chi_math::DynamicVector<double>(J + 1, 0.0);

  // ============================= Assembling system
  A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;
  for (size_t j = 1; j <= J; ++j)
  {
    A_[0][j] = lambdas_[j - 1];
    A_[j][j] = -lambdas_[j - 1];
    A_[j][0] = betas_[j - 1] / gen_time_;
  }

  q_.resize(J + 1, 0.0);
  q_[0] = source_strength_;

  // ============================= Initializing x
  // If there is a source and the reactivity is < 0 then
  // there exists a unique solution.
  if (source_strength_ > 0.0 and rho_ < 0.0)
  {
    const auto b_temp = -1.0 * q_;

    x_t_ = A_.Inverse() * b_temp;
  }
  // Otherwise we initialize the system as a critical system with
  // no source.
  else
  {
    auto A_temp = A_;
    auto b_temp = x_t_;
    b_temp.Set(0.0);
    b_temp[0] = 1.0;
    for (auto& val : A_temp[0])
      val = 0.0;
    A_temp[0][0] = 1.0;

    x_t_ = A_temp.Inverse() * b_temp;
  }

  chi::log.Log() << "Final: " << x_t_.PrintStr();
}

/**Execution function.*/
void TransientSolver::Execute() {}

/**Step function*/
void TransientSolver::Step()
{
  A_[0][0] = beta_ * (rho_ - 1.0) / gen_time_;

  const double theta = 1.0;
  const double inv_tau = theta * dt_;

  auto A_theta = I_ - inv_tau * A_;
  auto b_temp = x_t_ + inv_tau * q_;

  auto x_theta = A_theta.Inverse() * b_temp;

  x_tp1_ = x_t_ + (1.0 / theta) * (x_theta - x_t_);

  period_tph_ = dt_ / log(x_tp1_[0] / x_t_[0]);

  if (period_tph_ > 0.0 and period_tph_ > 1.0e6) period_tph_ = 1.0e6;
  if (period_tph_ < 0.0 and period_tph_ < -1.0e6) period_tph_ = -1.0e6;
}

/**Advance time values function.*/
void TransientSolver::Advance()
{
  time_ += dt_;
  x_t_ = x_tp1_;
}

} // namespace prk
