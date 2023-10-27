#include "QuadratureWarehouse.h"

namespace chi_math
{

QuadratureWarehouse& QuadratureWarehouse::GetInstance() noexcept
{
  static QuadratureWarehouse singleton;
  return singleton;
}

}