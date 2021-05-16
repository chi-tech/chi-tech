#ifndef LBS_STRUCTS_H
#define LBS_STRUCTS_H

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

  int max_iterations = 1000;
  double tolerance    = 1e-8;
  bool use_precursors = false;

  bool save_angular_flux = false;

  Options() = default;
};


/**Transport view of a cell*/
class CellLBSView
{
private:
  size_t phi_address;
  int num_nodes;
  int xs_mapping;
  std::vector<bool> face_local = {};
  int num_grps;
  int num_moms;

  std::vector<double> outflow;

public:
  CellLBSView(size_t in_phi_address,
              size_t in_num_nodes,
              size_t in_xs_mapping,
              const std::vector<bool>& face_local_flags,
              size_t num_G,
              size_t num_m,
              bool has_boundary_faces) :
    phi_address(in_phi_address),
    num_nodes(in_num_nodes),
    xs_mapping(in_xs_mapping),
    face_local(face_local_flags),
    num_grps(num_G),
    num_moms(num_m)
  {
    if (has_boundary_faces) outflow.assign(num_grps,0.0);
  }

  size_t MapDOF(int dof, int moment, int grp) const
  {
    return phi_address + dof * num_grps * num_moms + num_grps * moment + grp;
  }

  int NumNodes() const {return num_nodes;}
  int XSMapping() const {return xs_mapping;}
  bool IsFaceLocal(int f) const {return face_local[f];}
  void ZeroOutflow() {outflow.assign(num_grps,0.0);}
  void AddOutflow(int g, double value)
  {
    if (outflow.empty()) return;
    outflow[g] += value;
  }
};

}

#endif
