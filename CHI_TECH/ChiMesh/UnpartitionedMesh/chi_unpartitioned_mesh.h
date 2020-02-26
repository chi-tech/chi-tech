#ifndef chi_mesh_unpartitionedmesh_h
#define chi_mesh_unpartitionedmesh_h

#include "../chi_mesh.h"

#include <vtkCell.h>

//###################################################################
/**This object is intented for unpartitioned meshes that still require
 * partitioning.*/
class chi_mesh::UnpartitionedMesh
{
public:
  struct LightWeightFace
  {
    int neighbor=-1;
    std::vector<int> vertex_ids;
  };
  struct LightWeightCell
  {
    chi_mesh::Vertex centroid;
    std::vector<int> vertex_ids;
    std::vector<LightWeightFace> faces;
  };
private:
  std::vector<chi_mesh::Vertex*>  vertices;
  std::vector<LightWeightCell*>    raw_cells;

public:
  enum class ParallelMethod
  {
    ALL_FROM_HOME = 0,
    DIVIDE_WORK   = 1
  };
  struct Options
  {
    std::string file_name;
    ParallelMethod parallel_method = ParallelMethod::ALL_FROM_HOME;
  }mesh_options;

  LightWeightCell* CreateCellFromVTKPolyhedron(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKHexahedron(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKTetrahedron(vtkCell* vtk_cell);

  void ReadFromVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
};


#endif