#include "chi_surfaceremesher.h"

//################################################################### Default constructor
/**Default constructor.
\author Jan*/
CHI_SURFACEREMESHER::CHI_SURFACEREMESHER()
{
	this->baseSize = 0.01;
	this->absoluteMinumumSize = 0.005;
	this->removeInteriorFaces = true;
	this->keep2Dorientation = false;
	this->recalcNormal = false;
	this->precision = 0.0001;
}