#ifndef CHI_MESH_UNPARTITIONED_MESH_H
#define CHI_MESH_UNPARTITIONED_MESH_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

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
    std::vector<uint64_t> vertex_ids;
  };
  struct LightWeightCell
  {
    const chi_mesh::CellType type;
    chi_mesh::Vertex centroid;
    int material_id=-1;
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit
    LightWeightCell(chi_mesh::CellType in_type) : type(in_type) {}
  };

public:
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

  void BuildMeshConnectivity();
  void ComputeCentroidsAndCheckQuality();

  void ReadFromVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
  void ReadFromWavefrontOBJ(const Options& options);

  void ReadFromMsh(const Options& options);
};


#endif //CHI_MESH_UNPARTITIONED_MESH_H
