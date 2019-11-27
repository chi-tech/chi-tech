#include"chi_timer.h"


//############################################################################# Default constr
/** Default constructor.*/
ChiTimer::ChiTimer() noexcept
{
  startTime = std::chrono::steady_clock::now();
}
