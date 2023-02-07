#include "gs_context.h"

#include <petscksp.h>

namespace lbs
{

template<>
void GSContext<Mat, Vec>::ApplyInverseTransportOperator(int scope)
{

}

}//namespace lbs