#ifndef CHI_VOLUMEMESHER_H
#define CHI_VOLUMEMESHER_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

namespace chi_mesh
{
  enum VolumeMesherType
  {
    EXTRUDER      = 4,
    UNPARTITIONED = 6
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
}

//######################################################### Class def
/**Parent volume mesher class.*/
class chi_mesh::VolumeMesher
{
public:
  enum PartitionType
  {
//    KBA_STYLE_XY  = 1,
    KBA_STYLE_XYZ = 2,
    PARMETIS      = 3
  };
  struct VOLUME_MESHER_OPTIONS
  {
    bool         force_polygons = true;  //TODO: Remove this option
    bool         mesh_global    = false; //TODO: Remove this option
    int          partition_x    = 1;
    int          partition_y    = 1;
    int          partition_z    = 1;

    std::vector<double> xcuts;
    std::vector<double> ycuts;
    std::vector<double> zcuts;
    PartitionType partition_type = PARMETIS;
  };
  VOLUME_MESHER_OPTIONS options;
public:
  //01 Utils
  static
  void AddContinuumToRegion(MeshContinuumPtr& grid, Region& region);
  static
  void CreatePolygonCells(chi_mesh::SurfaceMesh* surface_mesh,
                          chi_mesh::MeshContinuumPtr& vol_continuum,
                          bool delete_surface_mesh_elements=false,
                          bool force_local=false);
  static
  void CreatePolygonCells(const chi_mesh::UnpartitionedMesh& umesh,
                          chi_mesh::MeshContinuumPtr& grid);
  static
  std::pair<int,int>  GetCellXYPartitionID(chi_mesh::Cell *cell);
  static
  std::tuple<int,int,int>
                      GetCellXYZPartitionID(chi_mesh::Cell *cell);
  static
  void                SetMatIDFromLogical(chi_mesh::LogicalVolume* log_vol,
                                          bool sense, int mat_id);
  static
  void                SetBndryIDFromLogical(chi_mesh::LogicalVolume* log_vol,
                                          bool sense, int bndry_id);
  static
  void                SetMatIDToAll(int mat_id);
  static
  void                SetupOrthogonalBoundaries();
  //02
  virtual void Execute();
  

};

#endif //CHI_VOLUMEMESHER_H