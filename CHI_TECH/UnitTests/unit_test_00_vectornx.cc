#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_tensorRNX.h"

#include "chi_log.h"

bool UnitTest_VectorNX()
{
  auto& log = ChiLog::GetInstance();


  {
    chi_math::Vector2 a(1.0);
    chi_math::Vector2 b(2.0);
    double s = 5.0;

    auto c=a;
    auto d=a;
    c*=b;
    d/=b;

    auto e=b;
    e.Normalize();

    auto f=a.Cross(b);

//    chi_math::VectorN<4> a4(1.0);
//    chi_math::VectorN<4> b4(2.0);
//
//    auto g = a4.Cross(b4);

    chi_math::Vector3 ihat(std::vector<double>{1.0,0.0});
    chi_math::Vector3 jhat(std::vector<double>{0.0,1.0,0.0,1.0});
    chi_math::Vector3 khat(std::vector<double>{0.0,0.0,1.0});

    auto g = ihat.Cross(jhat);

    log.Log(LOG_0) << "Vector a: " << a.PrintS();
    log.Log(LOG_0) << "Vector b: " << b.PrintS();
    log.Log(LOG_0) << "Scalar s: " << s;

    log.Log(LOG_0) << "a+b: " << (a+b).PrintS();
    log.Log(LOG_0) << "a-b: " << (a-b).PrintS();
    log.Log(LOG_0) << "a*b: " << (a*b).PrintS();
    log.Log(LOG_0) << "a/b: " << (a/b).PrintS();
    log.Log(LOG_0) << "a*=b: " << c.PrintS();
    log.Log(LOG_0) << "a/=b: " << d.PrintS();
    log.Log(LOG_0) << "s*a: " << (s*a).PrintS();
    log.Log(LOG_0) << "a*s: " << (a*s).PrintS();
    log.Log(LOG_0) << "b.Norm(): " << b.Norm();
    log.Log(LOG_0) << "b.NormSquare(): " << b.NormSquare();
    log.Log(LOG_0) << "b.Normalized(): " << b.Normalized().PrintS();
    log.Log(LOG_0) << "b.Normalize(): " << e.PrintS();
    log.Log(LOG_0) << "a.Cross(b): " << f.PrintS();
    log.Log(LOG_0) << "ihat: " << ihat.PrintS();
    log.Log(LOG_0) << "jhat: " << jhat.PrintS();
    log.Log(LOG_0) << "khat: " << khat.PrintS();
    log.Log(LOG_0) << "ihat.Cross(jhat): " << g.PrintS();
  }

  {
    chi_math::TensorRNX<2,5,double> t0(1.0);

    chi_math::VectorN<5> a(std::vector<double>{1.0,2.0,3.0,4.0,5.0});
    chi_math::VectorN<5> b(std::vector<double>{5.0,4.0,3.0,2.0,1.0});

    auto t1 = a.OTimes(b);

    auto c = a.Dot(t1);
    auto d = t1.Dot(a);

    log.Log(LOG_0) << "Vector a: " << a.PrintS();
    log.Log(LOG_0) << "Vector b: " << b.PrintS();

    log.Log(LOG_0) << "Tensor t0: \n" << t0.PrintS();
    log.Log(LOG_0) << "t1=a.OTimes(): \n" << t1.PrintS();
    log.Log(LOG_0) << "c=a.Dot(t1): \n" << c.PrintS();
    log.Log(LOG_0) << "d=t1.Dot(a): \n" << d.PrintS();
  }


  return true;
}