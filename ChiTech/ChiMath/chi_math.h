#ifndef CHI_MATH_H
#define CHI_MATH_H

#include "chi_math_incdef.h"

#include "Quadratures/quadrature.h"
#include "Quadratures/angular_quadrature_base.h"
#include "UnknownManager/unknown_manager.h"

#include <memory>


/**\defgroup LuaMath B Math*/

typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> MatDbl;

//######################################################### Class definition
/** This object handles the storage of math entities.

More information of solvers can be obtained from:
[HERE](http://www.maa.org/press/periodicals/loci/joma/iterative-methods-for-solving-iaxi-ibi-the-sor-method)
*/
class ChiMath
{
public:
	std::vector<chi_math::Quadrature*> quadratures;
	std::vector<std::shared_ptr<chi_math::AngularQuadrature>> angular_quadratures;

  static chi_math::UnknownManager UNITARY_UNKNOWN_MANAGER;
private:
  static ChiMath instance;
private:
	//00 Constructor
  ChiMath() noexcept
  {
    UNITARY_UNKNOWN_MANAGER.AddUnknown(chi_math::UnknownType::SCALAR);
  }

public:
  static ChiMath& GetInstance() noexcept
  {return instance;}
	//01 Utility

};

namespace chi_math
{
  /**Coordinate system type.*/
  enum class CoordinateSystemType
  {
    UNDEFINED   = 0,
    CARTESIAN   = 1,
  //CYLINDRICAL = 2,
  //SPHERICAL   = 3,
  };
  /**Spatial discretization type.*/
  enum class SpatialDiscretizationType
  {
    UNDEFINED                      = 0,
    FINITE_VOLUME                  = 1,
    PIECEWISE_LINEAR_CONTINUOUS    = 2,
    PIECEWISE_LINEAR_DISCONTINUOUS = 3,
    LAGRANGE_CONTINUOUS            = 4,
    LAGRANGE_DISCONTINUOUS         = 5
  };
  class SparseMatrix;

  class CDFSampler;
  int SampleCDF(double x, std::vector<double> cdf_bin);

  //02 Vector operations
  void   PrintVector(const VecDbl& x);
  void   Scale(VecDbl& x, const double& val);
  VecDbl VecMul(const VecDbl& x, const double& val);
  double Vec1Norm(const VecDbl& x);
  double Vec2Norm(const VecDbl& x);
  double VecInfinityNorm(const VecDbl& x);
  double VecPNorm(const VecDbl& x, const double& p);
  double Dot(const VecDbl& x, const VecDbl& y);

  //03 Matrix operations
  void   PrintMatrix(const MatDbl& A);
  void   Scale(MatDbl& A, const double& val);
  MatDbl Transpose(const MatDbl& A);
  void   SwapRow(size_t r1, size_t r2, MatDbl& A);
  void   SwapColumn(size_t c1, size_t c2, MatDbl& A);
  MatDbl MatMul(const MatDbl& A, const double c);
  VecDbl MatMul(const MatDbl& A, const VecDbl& x);
  MatDbl MatMul(const MatDbl& A, const MatDbl& B);
  MatDbl MatAdd(const MatDbl& A, const MatDbl& B);
  MatDbl MatSubtract(const MatDbl& A, const MatDbl& B);
  double Determinant(const MatDbl& A);
  MatDbl SubMatrix( const size_t r,
                    const size_t c,
                    const MatDbl& A );
  void   GaussElimination(MatDbl& A, VecDbl& b, int n);
  MatDbl InverseGEPivoting(const MatDbl& A);
  MatDbl Inverse(const MatDbl& A);

  double PowerIteration(const MatDbl& A,
                        VecDbl& e_vec, int max_it = 2000, double tol = 1.0e-13);

  //04 Unknown Managers
  class UnknownManager;
}


#endif
