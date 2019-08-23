#include"chi_surfaceremesher.h"

//############################################################################# Orient 2D
/**Determine the orientation of point a relative to point b and c.
 *
 * Synonomous with
\f[
	ab {\bullet} bc
\f]

but gives an indication of where $a$ lays relative to the line $bc$.

\return >0 if $a$ is right of $bc$, =0 if $abc$ are co-linear and <0 if $a$ is left of $bc$.

\image html Orient2D.png

\author Jan*/
double CHI_SURFACEREMESHER::Orient2D(float *a, float *b, float *c)
{
	Eigen::Vector3f A(a[0],a[1],a[2]);
	Eigen::Vector3f B(b[0],b[1],b[2]);
	Eigen::Vector3f C(c[0],c[1],c[2]);
	
	Eigen::Vector3f BA = B - A;
	Eigen::Vector3f CB = C - B;
	
	Eigen::Vector3f BAcrossCB = BA.cross(CB);
	
	return BAcrossCB(2);
}