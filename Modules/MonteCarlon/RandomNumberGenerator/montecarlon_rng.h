#ifndef _montecarlon_rng_h
#define _montecarlon_rng_h

#include "../chi_montecarlon.h"
#define R123_UNIFORM_FLOAT_STORE 1
#include<Random123/threefry.h>
#include<R123Uniform.hpp>

//#########################################################
/**Random number generator based on threefry.*/
class chi_montecarlon::RandomNumberGenerator
{
private:
  r123::Threefry2x64 generator;
  r123::Threefry2x64::ctr_type  ctr;
  r123::Threefry2x64::key_type key;

  bool    flipFlop;
  double  storedNumber;

public:
  //=================================== Default constructor
  /**Default constructor. Seed=0, Key=0*/
  RandomNumberGenerator()
  {
    ctr = {{0,0}};
    key = {{0,0}};

    flipFlop = false;
    storedNumber=0.5;
  }

  /**Constructor for a given seed and key.*/
  RandomNumberGenerator(uint64_t seed, uint64_t stream)
  {
    ctr = {{seed,seed}};
    key = {{stream,0}};

    flipFlop = false;
    storedNumber=0.5;
  }

  /**Returns a random number from the counter*/
  double Rand()
  {
    if (flipFlop==true)
    {
      flipFlop=false;
      return storedNumber;
    }
    else
    {
      flipFlop=true;
      ctr[0]++;
      ctr[1]++;
      r123::Threefry2x64::ctr_type rand = generator(ctr, key);
      storedNumber = r123::u01<double>(rand.v[0]);
      return         r123::u01<double>(rand.v[1]);
    }

  }
};


#endif