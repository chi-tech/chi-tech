#ifndef _surfacemesher_h
#define _surfacemesher_h

#include "../chi_mesh.h"

namespace chi_mesh
{
  enum class SurfaceMesherType
  {
    Passthrough = 1,
    Delaunay    = 2
  };
}

//###################################################################
/**Base class for surface meshers.*/
class chi_mesh::SurfaceMesher
{
public:
  const SurfaceMesherType type;
  int partitioning_x;
  int partitioning_y;
  bool export_loadbalance;

  std::vector<double> xcuts;
  std::vector<double> ycuts;
public:
  SurfaceMesher(SurfaceMesherType in_type) : type(in_type) {}

  virtual void Execute();

  void PrintLoadBalanceInfo();
};

#endif