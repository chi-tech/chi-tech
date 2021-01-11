#ifndef _chi_volumemesher_h
#define _chi_volumemesher_h

#include "../chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

//#define VOLUMEMESHER_LINEMESH1D 1
//#define VOLUMEMESHER_PREDEFINED2D 3
//#define VOLUMEMESHER_EXTRUDER 4

namespace chi_mesh
{
  enum VolumeMesherType
  {
    LINEMESH1D   = 1,
    PREDEFINED2D = 3,
    EXTRUDER     = 4,
    PREDEFINED3D = 5
  };
  enum VolumeMesherProperty
  {
    FORCE_POLYGONS      = 1,
    MESH_GLOBAL         = 2,
    PARTITION_Z         = 3,
    PARTITION_Y         = 4,
    PARTITION_X         = 5,
    CUTS_Z              = 6,
    CUTS_Y              = 7,
    CUTS_X              = 8,
    PARTITION_TYPE      = 9,
    EXTRUSION_LAYER     = 10,
    MATID_FROMLOGICAL   = 11,
    BNDRYID_FROMLOGICAL = 12
  };
};

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
public:
  std::vector<double> zcuts;

public:
  enum PartitionType
  {
    KBA_STYLE_XY  = 1,
    KBA_STYLE_XYZ = 2,
    PARMETIS      = 3
  };
  struct VOLUME_MESHER_OPTIONS
  {
    bool         force_polygons = true;
    bool         mesh_global    = false;
    int          partition_z    = 1;
    PartitionType partition_type = PARMETIS;
  };
  VOLUME_MESHER_OPTIONS options;
public:
  std::vector<chi_mesh::CellIndexMap*> cell_ordering;
  std::vector<chi_mesh::NodeIndexMap*> node_ordering;
  std::vector<int>                     reverse_node_ordering;
public:
  //01 Utils
  void AddContinuumToRegion(MeshContinuum* grid, Region& region);
  void CreatePolygonCells(chi_mesh::SurfaceMesh* surface_mesh,
                          chi_mesh::MeshContinuum* vol_continuum,
                          bool delete_surface_mesh_elements=false,
                          bool force_local=false);
  void GridFilterGhosts(chi_mesh::MeshContinuum *in_grid,
                        chi_mesh::MeshContinuum *out_grid);
  std::pair<int,int>  GetCellXYPartitionID(chi_mesh::Cell *cell);
  std::tuple<int,int,int>
                      GetCellXYZPartitionID(chi_mesh::Cell *cell);
  void                GetBoundaryCells(chi_mesh::MeshContinuum* vol_continuum);
  void                SetMatIDFromLogical(chi_mesh::LogicalVolume* log_vol,
                                          bool sense, int mat_id);
  void                SetBndryIDFromLogical(chi_mesh::LogicalVolume* log_vol,
                                          bool sense, int bndry_id);
  //02
  virtual void Execute();
  int          MapNode(int iref);
  int          ReverseMapNode(int i);
  

};

#endif