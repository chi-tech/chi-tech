#ifndef CHI_REMESHEDSURFACE_H
#define CHI_REMESHEDSURFACE_H
#include "../../../CHI_VECTOR/chi_vector.h"
#include"../../../CHI_GRAPHICS/CHI_SURFACE/chi_surface.h"
#include"../chi_surfaceremesher_structs.h"

using namespace chiSurfaceMeshing;
//############################################################################# Class Def
/** An object that is more advanced than a base CHI_SURFACE.

\author*/
class CHI_REMESHEDSURFACE : public CHI_SURFACE
{
public:
	CHI_VECTOR<CHI_EDGELOOP> edgeLoops;
	CHI_VECTOR<CHI_PATCH> patches;
	
public:
	CHI_REMESHEDSURFACE();
};

#endif