#ifndef CHI_TIMER_H
#define CHI_TIMER_H

#include <string>
#include <chrono>


//################################################################### CLASS DEF
namespace chi
{
  /** Timer object.*/
  class ChiTimer
  {
  private:
    std::chrono::steady_clock::time_point start_time_;

  public:
    //00
                ChiTimer() noexcept;
    //01
    void   		  Reset();
    double 		  GetTime() const;
    std::string GetTimeString() const;
    static std::string GetLocalDateTimeString();
  };
}

#endif

