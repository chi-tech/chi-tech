#include"chi_surfaceremesher.h"

//############################################################################# InCircle check
/** Checks if point d is within the circumcircle of triangle abc.
\author Jan*/
double CHI_SURFACEREMESHER::InCircle(float *a, float *b, float *c, float *d)
{
	Eigen::Matrix4f m;
	m(0,0) = a[0]; m(0,1) = a[1]; m(0,2) = a[0]*a[0] + a[1]*a[1]; m(0,3) = 1.0;
	m(1,0) = b[0]; m(1,1) = b[1]; m(1,2) = b[0]*b[0] + b[1]*b[1]; m(1,3) = 1.0;
	m(2,0) = c[0]; m(2,1) = c[1]; m(2,2) = c[0]*c[0] + c[1]*c[1]; m(2,3) = 1.0;
	m(3,0) = d[0]; m(3,1) = d[1]; m(3,2) = d[0]*d[0] + d[1]*d[1]; m(3,3) = 1.0;
	
	return m.determinant();
}