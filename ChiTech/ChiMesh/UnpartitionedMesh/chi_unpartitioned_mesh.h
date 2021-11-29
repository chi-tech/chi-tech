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
    std::vector<uint64_t> vertex_ids;
    bool has_neighbor = false;
    uint64_t neighbor=0;
  };
  struct LightWeightCell
  {
    const chi_mesh::CellType type;
    const chi_mesh::CellType sub_type;
    chi_mesh::Vertex centroid;
    int material_id=-1;
    std::vector<uint64_t> vertex_ids;
    std::vector<LightWeightFace> faces;

    explicit
    LightWeightCell(chi_mesh::CellType in_type,
                    chi_mesh::CellType in_sub_type) :
                    type(in_type),
                    sub_type(in_sub_type) {}
  };

public:
  std::vector<chi_mesh::Vertex>    vertices;
  std::vector<LightWeightCell*>    raw_cells;
  std::vector<LightWeightCell*>    raw_boundary_cells;
  std::vector<std::set<uint64_t>>  vertex_cell_subscriptions;

public:
  enum class ParallelMethod
  {
    ALL_FROM_HOME = 0,
    DIVIDE_WORK   = 1
  };
  struct Options
  {
    std::string file_name;
    std::string material_id_fieldname;
    std::string boundary_id_fieldname;
    double scale=1.0;
    ParallelMethod parallel_method = ParallelMethod::ALL_FROM_HOME;
  }mesh_options;

  struct BoundBox
  {
    double xmin=0.0, xmax=0.0,
           ymin=0.0, ymax=0.0,
           zmin=0.0, zmax=0.0;
  } bound_box;

  static LightWeightCell* CreateCellFromVTKPolyhedron(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKHexahedron(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKTetrahedron(vtkCell* vtk_cell);

  static LightWeightCell* CreateCellFromVTKPolygon(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKQuad(vtkCell* vtk_cell);
  static LightWeightCell* CreateCellFromVTKTriangle(vtkCell* vtk_cell);

  static LightWeightCell* CreateCellFromVTKLine(vtkCell* vtk_cell);

  static LightWeightCell* CreateCellFromVTKVertex(vtkCell* vtk_cell);

  void BuildMeshConnectivity();
  void ComputeCentroidsAndCheckQuality();

  void ReadFromVTU(const Options& options);
  void ReadFromEnsightGold(const Options& options);
  void ReadFromWavefrontOBJ(const Options& options);

  void ReadFromMsh(const Options& options);
};


#endif //CHI_MESH_UNPARTITIONED_MESH_H
