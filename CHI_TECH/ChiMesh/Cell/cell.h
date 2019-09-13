#ifndef _chi_cell_h
#define _chi_cell_h

#include"../chi_mesh.h"
#include <tuple>

//Appending cell types to namespace
namespace chi_mesh
{
  enum class CellType
  {
    GHOST = 0,
    SLAB = 1,
    SPHERICAL_SHELL = 2,
    CYLINDRICAL_ANNULUS = 3,
    TRIANGLE = 4,
    QUADRILATERAL = 5,
    POLYGON = 6,
    TETRAHEDRON = 7,
    HEXAHEDRON = 8,
    POLYHEDRON = 9
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

private:
  const CellType cell_type;
public:
  Cell(CellType in_cell_type) : cell_type(in_cell_type)
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

  virtual ~Cell() {}

public:
  virtual void FindBoundary2D(chi_mesh::Region* region) {}
  virtual bool CheckBoundary2D() {return true;}

  const CellType Type() {return cell_type;}
};


#endif