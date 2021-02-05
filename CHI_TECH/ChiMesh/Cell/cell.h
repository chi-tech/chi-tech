#ifndef CHI_CELL_H
#define CHI_CELL_H

#include"../chi_mesh.h"
#include <tuple>

//Appending cell types to namespace
namespace chi_mesh
{
enum class CellType
{
  GHOST = 0,
  SLAB = 1,
//  SPHERICAL_SHELL = 2,
//  CYLINDRICAL_ANNULUS = 3,
//  TRIANGLE = 4,
//  QUADRILATERAL = 5,
  POLYGON = 6,
//  TETRAHEDRON = 7,
//  HEXAHEDRON = 8,
  POLYHEDRON = 9,
};

//######################################################### Class def
/** In this paradigm a face is an object which largely
 * is considered to be planar (meaning all the vertices
 * lay in the same plane).*/
class CellFace
{
public:
  std::vector<uint64_t> vertex_ids; /// A list of the vertices
  Normal normal;                    ///< The average/geometric normal
  Vertex centroid;                  ///< The face centroid
  bool has_neighbor=false;          ///< Flag indicating whether face has a neighbor
  uint64_t neighbor_id=0;           ///< If face has neighbor, contains the global id. 0 otherwise.

public:
  bool IsNeighborLocal(chi_mesh::MeshContinuum& grid) const;
  int  GetNeighborPartitionID(chi_mesh::MeshContinuum& grid) const;
  int  GetNeighborLocalID(chi_mesh::MeshContinuum& grid) const;
  int  GetNeighborAssociatedFace(chi_mesh::MeshContinuum& grid) const;

public:
  double ComputeFaceArea(chi_mesh::MeshContinuum& grid) const;

};




//######################################################### Class def
/**Generic mesh cell object*/
class Cell
{
public:
  uint64_t global_id = 0;
  uint64_t local_id  = 0;
  uint64_t partition_id = 0;
  Vertex centroid;
  int material_id = -1;

  std::vector<uint64_t> vertex_ids;
  std::vector<CellFace> faces;

private:
  const CellType cell_type;

public:
  explicit Cell(CellType in_cell_type) : cell_type(in_cell_type) {}

  virtual ~Cell() = default;

public:
  CellType Type() const {return cell_type;}
};

}

#endif