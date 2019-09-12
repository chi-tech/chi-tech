#ifndef _chi_cell_h
#define _chi_cell_h

#include"../chi_mesh.h"
#include <tuple>

//######################################################### Class def
/**Generic mesh cell object*/
class chi_mesh::Cell
{
public:
  int cell_global_id;
  int cell_local_id;
  std::pair<int,int> xy_partition_indices;
  std::tuple<int,int,int> xyz_partition_indices;
  int partition_id;
  Vertex centroid;
  int material_id;
public:
  Cell()
  {
    cell_global_id = -1;
    cell_local_id = -1;
    xy_partition_indices.first  = 0;
    xy_partition_indices.second = 0;
    partition_id = -1;
    std::get<0>(xyz_partition_indices) = 0;
    std::get<1>(xyz_partition_indices) = 0;
    std::get<2>(xyz_partition_indices) = 0;

    material_id = -1;
  }

  virtual ~Cell()
  {}
public:
  virtual void FindBoundary2D(chi_mesh::Region* region)
  {}
  virtual bool CheckBoundary2D()
  {return true;}
};


#endif