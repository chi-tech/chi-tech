#ifndef _chi_domdecomp_h
#define _chi_domdecomp_h

#include <vector>
#include "../CHI_SURFACEMESHER/Triangle/triangle_mesher.h"

//###################################################################
/**Namespace that handles domain decomposition stuff.*/
namespace CHI_DOMDECOMP
{
  //01
  void Decompose2DDomain(int Px, int Py,
    chi_mesh::SurfaceMesherTriangle* mesher,
    chi_mesh::Region* region);
}


#endif