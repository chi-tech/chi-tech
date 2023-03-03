#ifndef CHI_LBS_ACCELERATION_H
#define CHI_LBS_ACCELERATION_H

#include <vector>

//################################################################### Fwd decls
namespace chi_physics
{
  class MultiGroupXSBase;
}


namespace lbs::acceleration
{

/**For both WGDSA and TGDSA use, we define this
* simplified data structure to hold multigroup diffusion coefficients and
* removal cross sections only for the relevant groups. E.g., for WGDSA it
* will only hold the cross sections for the groupset, and for TGDSA there will
* only be one group (All groups collapsed into 1).*/
struct Multigroup_D_and_sigR
{
  std::vector<double> Dg;
  std::vector<double> sigR;
};

enum class EnergyCollapseScheme
{
  JFULL    = 1, ///< Jacobi with full conv. of within-group scattering
  JPARTIAL = 2, ///< Jacobi with partially conv. of within-group scattering
};

struct TwoGridCollapsedInfo
{
  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum;
};

TwoGridCollapsedInfo
MakeTwoGridCollapsedInfo(const chi_physics::MultiGroupXSBase& xs,
                         EnergyCollapseScheme scheme);

}//namespace lbs::acceleration

#endif //CHI_LBS_ACCELERATION_H