#ifndef LBS_CURVILINEAR_SOLVER_H
#define LBS_CURVILINEAR_SOLVER_H

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"


namespace lbs_curvilinear
{
  class Solver;
}


/** A neutral particle transport solver in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class lbs_curvilinear::Solver : public lbs::SteadySolver
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
  Solver (const Solver&) = delete;
  Solver& operator= (const Solver&) = delete;

  /** Constructor. */
  Solver(const chi_math::CoordinateSystemType& coord_system_type,
         const std::string& in_text_name)
  : lbs::SteadySolver(in_text_name)
  , coord_system_type(coord_system_type)
  , discretization_secondary()
  {}

  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
private:
  std::shared_ptr<SweepChunk> SetSweepChunk(lbs::LBSGroupset& groupset) override;
};

#endif // LBS_CURVILINEAR_SOLVER_H
