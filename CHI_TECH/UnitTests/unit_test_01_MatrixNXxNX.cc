#include "ChiMath/chi_math_MatrixNXxNX.h"
#include "ChiMath/chi_math_vectorNX.h"
#include "chi_log.h"
#include <vector>
#include <iostream>

bool UnitTest_MatrixNXxNX()
{
  auto& log = ChiLog::GetInstance();
  {
      chi_math::MatrixNXxNX<3,double> a;
      chi_math::MatrixNXxNX<3,double> b(2.0);
      chi_math::MatrixNXxNX<3,double> c(3.0);
      double s = 5.0;
      std::vector<double> vec1 = {1,2,3};
      std::vector<double> vec2 = {4,5,6};
      std::vector<double> vec3 = {7,8,9};
      chi_math::VectorNX<3,double> vec4(7);

      log.Log(LOG_0) << "Matrix a: "<<"\n"<<a.PrintS();
      log.Log(LOG_0) << "Matrix b: "<<"\n"<<b.PrintS();
      log.Log(LOG_0) << "Matrix c: "<<"\n"<<c.PrintS();

      c = b;

      log.Log(LOG_0) << "assignment operator: " <<"\n"<< c.PrintS();
      log.Log(LOG_0) << "shift b+5: " <<"\n"<< b.Shift(5).PrintS();
      log.Log(LOG_0) << "b*vec4"<<"\n"<<(b*vec4).PrintS();
      log.Log(LOG_0) << "a+b: "<<"\n"<<(a+b).PrintS();
      log.Log(LOG_0) << "a+=b: " <<"\n"<<(a+=b).PrintS();
      log.Log(LOG_0) << "a-b: "<<"\n"<<(a-b).PrintS();
      log.Log(LOG_0) << "a-=b: "<<"\n"<<(a-=b).PrintS();
      log.Log(LOG_0) << "b*5: "<<"\n"<<(b*s).PrintS();
      log.Log(LOG_0) << "5*b: "<<"\n"<<(s*b).PrintS();
      log.Log(LOG_0) << "a*b: "<<"\n"<<(a*b).PrintS();
      log.Log(LOG_0) << "a*b component-wise" <<"\n"<< (a.MultComponentWise(b).PrintS());
      log.Log(LOG_0) << "b/5" <<"\n"<<(b/s).PrintS();
      log.Log(LOG_0) << "a/b" <<"\n"<<(a.DivComponentWise(b)).PrintS();
      log.Log(LOG_0) << "get (1,1) of a" <<"\n"<<a.GetIJ(1,1);
      
      a.SetIJ(1,1,5.0);
      log.Log(LOG_0) << "setIJ (1,1) of a with 5" <<"\n" <<a.PrintS();
      a.AddIJ(1,1,5.0);
      log.Log(LOG_0) << "addIJ (1,1) of a with 5" <<"\n" <<a.PrintS();
      a.SetDiagonalVec(vec1);
      log.Log(LOG_0) << "diagonal vec using {1,2,3}" <<"\n" <<a.PrintS();
      a.SetRowIVec(0,vec1);
      log.Log(LOG_0) << "row vec using {1,2,3}" <<"\n" <<a.PrintS();
      a.SetColJVec(0,vec1);
      log.Log(LOG_0) << "column vec using {1,2,3}" <<"\n" <<a.PrintS();
      
      a.SetRowIVec(0,{1,2,3});
      a.SetRowIVec(1,{4,5,6});
      a.SetRowIVec(2,{7,8,4});

      log.Log(LOG_0) << "finding determinant of matrix"<<a.Det();
      log.Log(LOG_0) << "finding transpose"<<"\n"<<a.Transpose().PrintS();
      log.Log(LOG_0) << "finding inverse"<<"\n"<<a.Inverse().PrintS();
      log.Log(LOG_0) << "norm of a: "<<a.Norm();
  }
  return true;
}



