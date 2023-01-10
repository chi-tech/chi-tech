#include "ChiMath/chi_math_vectorNX.h"
#include "ChiMath/chi_math_tensorRNX.h"

#include <iostream>

int main()
{
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

    chi_math::Vector3 ihat(std::vector<double>{1.0,0.0});
    chi_math::Vector3 jhat(std::vector<double>{0.0,1.0,0.0,1.0});
    chi_math::Vector3 khat(std::vector<double>{0.0,0.0,1.0});

    auto g = ihat.Cross(jhat);

    std::cout << "Vector a: " << a.PrintS() << "\n";
    std::cout << "Vector b: " << b.PrintS() << "\n";
    std::cout << "Scalar s: " << s << "\n";

    std::cout << "a+b: " << (a+b).PrintS() << "\n";
    std::cout << "a-b: " << (a-b).PrintS() << "\n";
    std::cout << "a*b: " << (a*b).PrintS() << "\n";
    std::cout << "a/b: " << (a/b).PrintS() << "\n";
    std::cout << "a*=b: " << c.PrintS() << "\n";
    std::cout << "a/=b: " << d.PrintS() << "\n";
    std::cout << "s*a: " << (s*a).PrintS() << "\n";
    std::cout << "a*s: " << (a*s).PrintS() << "\n";
    std::cout << "b.Norm(): " << b.Norm() << "\n";
    std::cout << "b.NormSquare(): " << b.NormSquare() << "\n";
    std::cout << "b.Normalized(): " << b.Normalized().PrintS() << "\n";
    std::cout << "b.Normalize(): " << e.PrintS() << "\n";
    std::cout << "a.Cross(b): " << f.PrintS() << "\n";
    std::cout << "ihat: " << ihat.PrintS() << "\n";
    std::cout << "jhat: " << jhat.PrintS() << "\n";
    std::cout << "khat: " << khat.PrintS() << "\n";
    std::cout << "ihat.Cross(jhat): " << g.PrintS() << "\n";
  }

  {
    chi_math::TensorRNX<2,5,double> t0(1.0);

    chi_math::VectorN<5> a(std::vector<double>{1.0,2.0,3.0,4.0,5.0});
    chi_math::VectorN<5> b(std::vector<double>{5.0,4.0,3.0,2.0,1.0});

    auto t1 = a.OTimes(b);

    auto c = a.Dot(t1);
    auto d = t1.Dot(a);

    std::cout << "Vector a: " << a.PrintS() << "\n";
    std::cout << "Vector b: " << b.PrintS() << "\n";

    std::cout << "Tensor t0: \n" << t0.PrintS() << "\n";
    std::cout << "t1=a.OTimes(): \n" << t1.PrintS() << "\n";
    std::cout << "c=a.Dot(t1): \n" << c.PrintS() << "\n";
    std::cout << "d=t1.Dot(a): \n" << d.PrintS() << "\n";
  }
  return 0;
}