#include "chi_math.h"
#include <assert.h>

//######################################################### Print Vector
/** Prints the Vector.*/
void ChiMath::PrintVector(const VecDbl &x)
{
  for (auto& xi : x)
    std::cout << xi << ' ';
  std::cout << std::endl;
}

//#########################################################  Scale
/** Scales a vector in place by constant.*/
void ChiMath::Scale(VecDbl &x, const double &val)
{
  for (double& xi : x)
    xi *= val;
}

//######################################################### Dot Product
/** Computes the dot product of two vectors.
 *
 * \f[
 * \mathrm{a} \cdot \mathrm{b}=\sum_{i=1}^{n} a_{i} b_{i}
 * \f]
 * */
double ChiMath::Dot(const VecDbl &x, const VecDbl &y)
{
  // Error Checks
  assert(x.size() > 0);
  assert(y.size() > 0);
  assert(x.size() == y.size());
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for(size_t i = 0; i != n; i++)
    val += x[i] * y[i];

  return val;
}

//######################################################### Multiply with const
/** Multiplies the vector with a constant and returns result.*/
VecDbl ChiMath::VecMul(const VecDbl &x, const double &val)
{
  size_t n = x.size();
  VecDbl y(n);

  for (size_t i = 0; i != n; ++i)
    y[i] = val * x[i];

  return y;
}

//######################################################### Norms
/** Returns the 1-norm. Also known as the Taxicab or Manhattan norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{1}=\sum_{i=1}^{n}\left|x_{i}\right|
 * \f]
 *
 * */
 double ChiMath::Vec1Norm(const VecDbl &x)
 {
   // Local Variables
   size_t n = x.size();
   double val = 0.;

   for(size_t i = 0; i != n; i++)
     val += std::fabs(x[i]);

   return val;
 }


/** Returns the 2-norm. Also known as the Euclidian or Frobenius norm.
 *
 * \f[
 * \|\boldsymbol{x}\|_{2}=\sqrt{x_{1}^{2}+\cdots+x_{n}^{2}}
 * \f]
 *
 * */
double ChiMath::Vec2Norm(const VecDbl &x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.;

  for(size_t i = 0; i != n; i++)
    val += x[i]*x[i];

  return std::sqrt(val);
}

/** Returns the infinity-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{\infty}=\max \left(\left|x_{1}\right|, \ldots,\left|x_{n}\right|\right)
 * \f]
 *
 * */
double ChiMath::VecInfinityNorm(const VecDbl &x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for(size_t i = 0; i != n; i++)
    val += std::max(std::fabs(x[i]), val);

  return val;
}

/** Returns the p-norm.
 *
 * \f[
 * \|\mathbf{x}\|_{p}=\left(\sum_{i=1}^{n}\left|x_{i}\right|^{p}\right)^{1 / p}
 * \f]
 *
 * */
double ChiMath::VecPNorm(const VecDbl &x, const double& p)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.;

  for(size_t i = 0; i != n; i++)
    val += std::pow(std::fabs(x[i]), p);

  return std::pow(val, 1./p);
}

