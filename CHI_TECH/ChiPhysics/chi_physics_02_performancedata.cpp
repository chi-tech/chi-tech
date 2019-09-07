#include "chi_physics.h"
#include<iostream>
#include<fstream>

//############################################################################# Print performance data
/** Prints the performance data array to file

\author Jan*/
void ChiPhysics::PrintPerformanceData(char* fileName)
{
	  std::ofstream file;
	  file.open(fileName);
//	  file.setf(std::ios::scientific);
//	  file.precision(12);
//
//	  //======================================================= Write ofstream
//	  for (long int i=0;i<10000;i++)
//	  {
//		  file << this->performanceData[i] << std::endl;
//	  }



	  file.close();


}
