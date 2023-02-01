#ifndef LBS_CURVILINEAR_SOLVER_H
#define LBS_CURVILINEAR_SOLVER_H

#include "LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"


namespace lbs_curvilinear
{

/** A neutral particle transport solver in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class SteadyStateSolver : public lbs::SteadyStateSolver
{
//  Attributes
private:
  /** Coordinate system type. */
  chi_math::CoordinateSystemType coord_system_type;
  /** Discretisation pointer to matrices of the secondary cell view
   *  (matrices of the primary cell view forwarded to the base class). */
  std::shared_ptr<chi_math::SpatialDiscretization> discretization_secondary;

//  Methods
public:
  SteadyStateSolver (const SteadyStateSolver&) = delete;
  SteadyStateSolver& operator= (const SteadyStateSolver&) = delete;

  /** Constructor. */
  SteadyStateSolver(const chi_math::CoordinateSystemType& coord_system_type,
                    const std::string& in_text_name)
  : lbs::SteadyStateSolver(in_text_name)
  , coord_system_type(coord_system_type)
  , discretization_secondary()
  {}

  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
private:
  std::shared_ptr<SweepChunk> SetSweepChunk(lbs::LBSGroupset& groupset) override;
};

}//namespace lbs_curvilinear

#endif // LBS_CURVILINEAR_SOLVER_H
