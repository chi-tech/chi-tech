#ifndef chi_mesh_unpartitionedmesh_h
#define chi_mesh_unpartitionedmesh_h

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

#include <vtkCell.h>

//###################################################################
/**This object is intented for unpartitioned meshes that still require
 * partitioning.*/
class chi_mesh::UnpartitionedMesh
{
  friend class VolumeMesherPredefined3D;
  friend class VolumeMesherPredefinedUnpartitioned;
  friend void CreateUnpartitioned1DOrthoMesh(std::vector<double>& vertices);
  friend void CreateUnpartitioned2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                             std::vector<double>& vertices_1d_y);
  friend void CreateUnpartitioned3DOrthoMesh(std::vector<double>& vertices_1d_x,
                                             std::vector<double>& vertices_1d_y,
                                             std::vector<double>& vertices_1d_z);
private:
  struct LightWeightFace
  {
    int neighbor=-1;
    std::vector<int> vertex_ids;
  };
  struct LightWeightCell
  {
    const chi_mesh::CellType type;
    chi_mesh::Vertex centroid;
    int material_id=-1;
    std::vector<int> vertex_ids;
    std::vector<LightWeightFace> faces;

    LightWeightCell(chi_mesh::CellType in_type) : type(in_type) {}
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
    double scale=1.0;
    ParallelMethod parallel_method = ParallelMethod::ALL_FROM_HOME;
  }mesh_options;

  struct BoundBox
  {
    double xmin=0.0, xmax=0.0,
           ymin=0.0, ymax=0.0,
           zmin=0.0, zmax=0.0;
  } bound_box;

  LightWeightCell* CreateCellFromVTKPolyhedron(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKHexahedron(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKTetrahedron(vtkCell* vtk_cell);

  LightWeightCell* CreateCellFromVTKPolygon(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKQuad(vtkCell* vtk_cell);
  LightWeightCell* CreateCellFromVTKTriangle(vtkCell* vtk_cell);

  void ReadFromVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
};


#endif