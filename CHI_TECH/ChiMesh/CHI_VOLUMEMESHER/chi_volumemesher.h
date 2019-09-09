#ifndef _chi_volumemesher_h
#define _chi_volumemesher_h

#include "../chi_mesh.h"
#include "../CHI_CELL/cell.h"

#define VOLUMEMESHER_LINEMESH1D 1
#define VOLUMEMESHER_PREDEFINED2D 3
#define VOLUMEMESHER_EXTRUDER 4

struct chi_mesh::CellIndexMap
{
  int mapped_from;
  int mapped_to;
  int mapped_level;

  CellIndexMap()
  {
    mapped_from  = -1;
    mapped_to    = -1;
    mapped_level = -1;
  }
  CellIndexMap(int from,int to)
  {
    mapped_from  = from;
    mapped_to    = to;
    mapped_level = -1;
  }
};

//######################################################### Class def
/**Parent volume mesher class.*/
class chi_mesh::VolumeMesher
{
private:
  std::vector<double> zcuts;

public:
  struct VOLUME_MESHER_OPTIONS
  {
    bool force_polygons;
    bool mesh_global;
    int  partition_z;

    VOLUME_MESHER_OPTIONS()
    {
      force_polygons = true;
      mesh_global = false;
      partition_z = 1;
    }
  };
  VOLUME_MESHER_OPTIONS options;
public:
  std::vector<chi_mesh::CellIndexMap*> cell_ordering;
  std::vector<chi_mesh::NodeIndexMap*> node_ordering;
  std::vector<int>                     reverse_node_ordering;
public:
  //01 Utils
  void                CreateTriangleCells(chi_mesh::SurfaceMesh* surface_mesh,
                           chi_mesh::MeshContinuum* vol_continuum);
  void                CreatePolygonCells(chi_mesh::SurfaceMesh* surface_mesh,
                          chi_mesh::MeshContinuum* vol_continuum);
  std::pair<int,int>  GetCellXYPartitionID(chi_mesh::Cell *cell);
  std::tuple<int,int,int>
                      GetCellXYZPartitionID(chi_mesh::Cell *cell);
  void                GetBoundaryCells(chi_mesh::MeshContinuum* vol_continuum);
  void                SetMatIDFromLogical(chi_mesh::LogicalVolume* log_vol,
                                          bool sense, int mat_id);
  //02
  virtual void Execute();
  int          MapNode(int iref);
  int          ReverseMapNode(int i);

  //03
  void         ReOrderDOF(chi_mesh::MeshContinuum* vol_continuum);
  void         ReOrderDOF_CuthillMckee(chi_mesh::MeshContinuum* vol_continuum);
  bool         IsInList(std::vector<chi_mesh::NodeIndexMap*>* list, int index);

  int          MapLevel(std::vector<chi_mesh::NodeIndexMap*>* used_list,
                        std::vector<chi_mesh::NodeIndexMap*>* unused_cells,
                        chi_mesh::MeshContinuum* vol_continuum,
                        int seeded_index,
                        int level);


};

#endif