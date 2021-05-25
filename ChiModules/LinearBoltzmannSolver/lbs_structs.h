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
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL   = 3,
  TWOD_CARTESIAN   = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6
};

/**Struct for storing LBS options.*/
struct Options
{
  typedef chi_math::SpatialDiscretizationType SDMType;

  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  SDMType sd_type = SDMType::UNDEFINED;
  unsigned int scattering_order=1;
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
  int num_grps;
  int num_grps_moms;
  int xs_mapping;
  std::vector<bool> face_local_flags = {};
  std::vector<double> outflow;

public:
  CellLBSView(size_t in_phi_address,
              int in_num_nodes,
              int in_num_grps,
              int in_num_moms,
              int in_xs_mapping,
              const std::vector<bool>& in_face_local_flags,
              bool cell_on_boundary) :
    phi_address(in_phi_address),
    num_nodes(in_num_nodes),
    num_grps(in_num_grps),
    num_grps_moms(in_num_grps*in_num_moms),
    xs_mapping(0),
    face_local_flags(in_face_local_flags)
  {
    if (cell_on_boundary)
      outflow.resize(num_grps,0.0);
  }

  size_t MapDOF(int node, int moment, int grp) const
  {
    return phi_address + node * num_grps_moms + num_grps * moment + grp;
  }

  int XSMapping() const {return xs_mapping;}

  bool IsFaceLocal(int f) const {return face_local_flags[f];}

  int NumNodes() const {return num_nodes;}

  void ZeroOutflow() {outflow.assign(outflow.size(),0.0);}
  void AddOutflow(int g, double intS_mu_psi) {outflow[g] += intS_mu_psi;}
  double GetOutflow(int g) const
  {
    if (g<outflow.size()) return outflow[g];
    else return 0.0;
  }
};

}

#endif
