#include "quadrature.h"

#include <cassert>


void
chi_math::Quadrature::
  SetRange(const std::pair<double, double>& in_range)
{
  const auto& old_range = range;
  const auto& new_range = in_range;

  const double h_new = new_range.second - new_range.first;
  const double h_old = old_range.second - old_range.first;

  if (h_new <= 0.0 or h_old <= 0.0)
    throw std::invalid_argument("Quadrature::"+std::string(__FUNCTION__)+
                                ": called with negative or zero ranges.");

  if (qpoints.empty())
    throw std::invalid_argument("Quadrature::"+std::string(__FUNCTION__)+
                                ": called with no abscissae initialized.");

  const double scale_factor = h_new/h_old;

  for (unsigned int i=0; i < qpoints.size(); ++i)
  {
    qpoints[i](0) = new_range.first +
                    (qpoints[i][0] - old_range.first) * scale_factor;

    weights[i] *= scale_factor;
  }

  range = in_range;
}
