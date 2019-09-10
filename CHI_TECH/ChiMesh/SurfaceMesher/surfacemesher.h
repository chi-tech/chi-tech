#ifndef _surfacemesher_h
#define _surfacemesher_h

#include "../chi_mesh.h"
class chi_mesh::SurfaceMesher
{
public:
  int partitioning_x;
  int partitioning_y;
  bool export_loadbalance;

  std::vector<double> xcuts;
  std::vector<double> ycuts;
public:
  virtual void Execute();

  void PrintLoadBalanceInfo();
};

#endif