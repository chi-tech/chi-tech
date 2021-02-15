#ifndef LBS_STRUCTS_H
#define LBS_STRUCTS_H

#define PARTITION_METHOD_SERIAL        1
#define PARTITION_METHOD_FROM_SURFACE  2

#include "ChiMath/chi_math.h"

namespace LinearBoltzmann
{

enum class GeometryType
{
  NO_GEOMETRY_SET  = 0,
  ONED_SLAB        = 1,
  ONED_SPHERICAL   = 2,
  TWOD_CARTESIAN   = 3,
  TWOD_CYLINDRICAL = 4,
  THREED_CARTESIAN = 5
};

/**Struct for storing LBS options.*/
struct Options
{
  typedef chi_math::SpatialDiscretizationType SDMType;

  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  SDMType sd_type = SDMType::UNDEFINED;
  int  scattering_order=1;
  int  sweep_eager_limit= 32000;;

  bool read_restart_data=false;
  std::string read_restart_folder_name = std::string("YRestart");
  std::string read_restart_file_base   = std::string("restart");

  bool write_restart_data=false;
  std::string write_restart_folder_name = std::string("YRestart");
  std::string write_restart_file_base   = std::string("restart");
  double write_restart_interval = 30.0;

  Options() = default;
};


/**Transport view of a cell*/
class CellLBSView
{
public:
  int dof_phi_map_start = 0;
  int dofs = 0;
  int xs_id = 0;
  std::vector<bool> face_local = {};

private:
  int num_grps = 0;
  int num_moms = 0;

public:
  CellLBSView(int in_dofs, int num_G, int num_m)
  {
    dof_phi_map_start = -1;
    dofs = in_dofs;
    num_grps = num_G;
    num_moms = num_m;
  }

  int MapDOF(int dof, int moment, int grp) const
  {
    return dof_phi_map_start + dof*num_grps*num_moms + num_grps*moment + grp;
  }
};

}

#endif
