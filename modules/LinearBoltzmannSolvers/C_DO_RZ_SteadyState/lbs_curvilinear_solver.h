#ifndef LBS_CURVILINEAR_SOLVER_H
#define LBS_CURVILINEAR_SOLVER_H

#include "B_DO_SteadyState/lbs_DO_steady_state.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"


namespace lbs_curvilinear
{

/** A neutral particle transport solver in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class DiscOrdSteadyStateSolver : public lbs::DiscOrdSteadyStateSolver
{
//  Attributes
private:
  /** Coordinate system type. */
  chi_math::CoordinateSystemType coord_system_type_;
  /** Discretisation pointer to matrices of the secondary cell view
   *  (matrices of the primary cell view forwarded to the base class). */
  std::shared_ptr<chi_math::SpatialDiscretization> discretization_secondary_;

//  Methods
public:
  DiscOrdSteadyStateSolver (const DiscOrdSteadyStateSolver&) = delete;
  DiscOrdSteadyStateSolver& operator= (const DiscOrdSteadyStateSolver&) = delete;

  /** Constructor. */
  DiscOrdSteadyStateSolver(const chi_math::CoordinateSystemType& coord_system_type,
                           const std::string& in_text_name)
  : lbs::DiscOrdSteadyStateSolver(in_text_name)
  , coord_system_type_(coord_system_type)
  , discretization_secondary_()
  {}

  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
private:
  std::shared_ptr<SweepChunk> SetSweepChunk(lbs::LBSGroupset& groupset) override;
};

}//namespace lbs_curvilinear

#endif // LBS_CURVILINEAR_SOLVER_H
