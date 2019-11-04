#ifndef _cell_newbase_h
#define _cell_newbase_h

#include "ChiMesh/chi_mesh.h"
#include "cell.h"
namespace chi_mesh
{
//######################################################### Class def
/** In this paradigm a face is an object which largely
 * is considered to be planar (meaning all the vertices
 * lay in the same plane).*/
class CellFace
{
public:
  std::vector<int> vertex_ids;   /// A list of the vertices
  Normal normal;                 /// The average/geometric normal
  Vertex centroid;               /// The face centroid
  int neighbor;                  /// Neigboring cell index (<0 indicates bndry)

  CellFace()
  {
    neighbor = -1;
  }
};

//######################################################### Class def
/**Generic mesh cell object*/
class CellBase : public Cell
{
public:
//  int cell_global_id;
//  int cell_local_id;
//  std::pair<int,int> xy_partition_indices;
//  std::tuple<int,int,int> xyz_partition_indices;
//  int partition_id;
//  Vertex centroid;
//  int material_id;

private:
//  const CellType cell_type;
    const CellType cell_type2;
public:
//  explicit Cell(CellType in_cell_type) : cell_type(in_cell_type)
//  {
//    cell_global_id = -1;
//    cell_local_id = -1;
//    xy_partition_indices.first  = 0;
//    xy_partition_indices.second = 0;
//    partition_id = -1;
//    std::get<0>(xyz_partition_indices) = 0;
//    std::get<1>(xyz_partition_indices) = 0;
//    std::get<2>(xyz_partition_indices) = 0;
//
//    material_id = -1;
//  }
//
//  virtual ~CellBase() {}
//
//public:
//  virtual void FindBoundary2D(chi_mesh::Region* region) {}
//  virtual bool CheckBoundary2D() {return true;}
//
  const CellType Type2() {return cell_type2;}

  std::vector<int> vertex_ids;
  std::vector<CellFace> faces;

  CellBase(CellType in_cell_type) :
    Cell(CellType::CELL_NEWBASE),
    cell_type2(in_cell_type)
    {}
};

}

#endif