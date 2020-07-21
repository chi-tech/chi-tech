#ifndef _lbs_structs_h
#define _lbs_structs_h

#define PARTITION_METHOD_SERIAL        1
#define PARTITION_METHOD_FROM_SURFACE  2

namespace LinearBoltzman
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

/**Struct for storing NPT options.*/
struct Options
{
  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  int  scattering_order;
  int  partition_method;
  int  sweep_eager_limit;

  bool read_restart_data;
  std::string read_restart_folder_name;
  std::string read_restart_file_base;

  bool write_restart_data;
  std::string write_restart_folder_name;
  std::string write_restart_file_base;
  double write_restart_interval;

  Options()
  {
    scattering_order = 0;
    partition_method = PARTITION_METHOD_SERIAL;
    sweep_eager_limit= 32000;

    read_restart_data = false;
    read_restart_folder_name = std::string("YRestart");
    read_restart_file_base   = std::string("restart");

    write_restart_data = false;
    write_restart_folder_name = std::string("YRestart");
    write_restart_file_base   = std::string("restart");
    write_restart_interval = 30.0;
  }
};

class CellViewBase
{
public:
  std::vector<int> node_dof_mapping;
};



/**Transport view of a cell*/
class CellViewFull : public CellViewBase
{
public:
  int dof_phi_map_start;
  int dofs;
  int xs_id;
  std::vector<bool> face_local;

private:
  int num_grps;
  int num_moms;

public:
  CellViewFull(int in_dofs, int num_G, int num_m)
  {
    dof_phi_map_start = -1;
    dofs = in_dofs;
    num_grps = num_G;
    num_moms = num_m;
  }

  int MapDOF(int dof, int moment, int grp)
  {
    return dof_phi_map_start + dof*num_grps*num_moms + num_grps*moment + grp;
  }
};

}

#endif
