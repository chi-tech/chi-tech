#include "delaunay_mesher.h"

/**Default constructor.*/
chi_mesh::SurfaceMesherDelaunay::SurfaceMesherDelaunay() :
  SurfaceMesher(SurfaceMesherType::Delaunay)
{
  tolerance = 1.0e-5;
}