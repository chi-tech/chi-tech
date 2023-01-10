#ifndef _surfacemesher_h
#define _surfacemesher_h

#include "../chi_mesh.h"

namespace chi_mesh
{
  enum class SurfaceMesherType
  {
    Predefined = 1
  };
}

//###################################################################
/**Base class for surface meshers.*/
class chi_mesh::SurfaceMesher
{
public:
  const SurfaceMesherType type;

  std::vector<double> xcuts;
  std::vector<double> ycuts;
public:
  explicit SurfaceMesher(SurfaceMesherType in_type) : type(in_type) {}

  virtual void Execute();

};

#endif