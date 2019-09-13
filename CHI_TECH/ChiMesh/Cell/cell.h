#ifndef _chi_cell_h
#define _chi_cell_h

#include"../chi_mesh.h"
#include <tuple>

namespace chi_mesh
{
  enum CellTypes
  {
    GHOST_CELL = 0,
    SLAB_CELL = 1,
    SPHERICAL_SHELL_CELL = 2,
    CYLINDRICAL_ANNULUS_CELL = 3,
    TRIANGLE_CELL = 4,
    QUADRILATERAL_CELL = 5,
    POLYGON_CELL = 6,
    TETRAHEDRON_CELL = 7,
    HEXAHEDRON_CELL = 8,
    POLYHEDRON_CELL = 9
  };
}



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

protected:
  CellTypes cell_type;
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
  {
    cell_type = GHOST_CELL;
  }
public:
  virtual void FindBoundary2D(chi_mesh::Region* region)
  {}
  virtual bool CheckBoundary2D()
  {return true;}
  CellTypes Type()
  {
    return cell_type;
  }
};


#endif