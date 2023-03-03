#ifndef LBS_STRUCTS_H
#define LBS_STRUCTS_H

#include "ChiMath/chi_math.h"
#include "ChiPhysics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include <functional>

namespace lbs
{

typedef std::vector<size_t> DirIDs; ///< Direction-IDs
typedef std::vector<DirIDs> UniqueSOGroupings;
typedef std::map<size_t, size_t> DirIDToSOMap;

typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> MatDbl;
typedef std::vector<chi_mesh::Vector3> VecVec3;
typedef std::vector<VecVec3> MatVec3;

typedef std::shared_ptr<chi_physics::MultiGroupXSBase> XSPtr;
typedef std::shared_ptr<chi_physics::IsotropicMultiGrpSource> IsotropicSrcPtr;

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

enum class AngleAggregationType
{
  UNDEFINED = 0,
  SINGLE = 1,
  POLAR = 2,
  AZIMUTHAL = 3,
};

enum class BoundaryType
{
  VACUUM = 1,
  INCIDENT_ISOTROPIC = 2,
  REFLECTING = 3
};

struct BoundaryPreference
{
  BoundaryType type;
  std::vector<double> isotropic_mg_source;
  std::string source_function;
};

enum SourceFlags : int
{
  NO_FLAGS_SET              = 0,
  APPLY_FIXED_SOURCES       = (1 << 0),
  APPLY_WGS_SCATTER_SOURCES = (1 << 1),
  APPLY_AGS_SCATTER_SOURCES = (1 << 2),
  APPLY_WGS_FISSION_SOURCES = (1 << 3),
  APPLY_AGS_FISSION_SOURCES = (1 << 4),
  SUPPRESS_WG_SCATTER       = (1 << 5)
};

inline SourceFlags operator|(const SourceFlags f1,
                             const SourceFlags f2)
{
  return static_cast<SourceFlags>(static_cast<int>(f1) |
                                  static_cast<int>(f2));
}

enum class PhiSTLOption
{
  PHI_OLD = 1,
  PHI_NEW = 2
};

class LBSGroupset;
typedef std::function<void(LBSGroupset&               groupset,
                           std::vector<double>&       destination_q,
                           const std::vector<double>& phi,
                           SourceFlags                source_flags)>
                      SetSourceFunction;

class AGSSchemeEntry;

/**Struct for storing LBS options.*/
struct Options
{
  typedef chi_math::SpatialDiscretizationType SDMType;

  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  SDMType sd_type = SDMType::UNDEFINED;
  unsigned int scattering_order=1;
  int  sweep_eager_limit= 32000; //see chiLBSSetProperty documentation

  bool read_restart_data=false;
  std::string read_restart_folder_name = std::string("YRestart");
  std::string read_restart_file_base   = std::string("restart");

  bool write_restart_data=false;
  std::string write_restart_folder_name = std::string("YRestart");
  std::string write_restart_file_base   = std::string("restart");
  double write_restart_interval = 30.0;

  bool use_precursors = false;
  bool use_src_moments = false;

  bool save_angular_flux = false;

  bool verbose_inner_iterations = true;
  bool verbose_ags_iterations = false;
  bool verbose_outer_iterations = true;

  Options() = default;

  std::vector<AGSSchemeEntry> ags_scheme;
};


/**Transport view of a cell*/
class CellLBSView
{
private:
  size_t phi_address;
  int num_nodes;
  int num_grps;
  int num_grps_moms;
  const chi_physics::MultiGroupXSBase* xs;
  double volume;
  std::vector<bool> face_local_flags = {};
  std::vector<double> outflow;

public:
  CellLBSView(size_t in_phi_address,
              int in_num_nodes,
              int in_num_grps,
              int in_num_moms,
              const chi_physics::MultiGroupXSBase& in_xs_mapping,
              double in_volume,
              const std::vector<bool>& in_face_local_flags,
              bool cell_on_boundary) :
    phi_address(in_phi_address),
    num_nodes(in_num_nodes),
    num_grps(in_num_grps),
    num_grps_moms(in_num_grps*in_num_moms),
    xs(&in_xs_mapping),
    volume(in_volume),
    face_local_flags(in_face_local_flags)
  {
    if (cell_on_boundary)
      outflow.resize(num_grps,0.0);
  }

  size_t MapDOF(int node, int moment, int grp) const
  {
    return phi_address + node * num_grps_moms + num_grps * moment + grp;
  }

  const chi_physics::MultiGroupXSBase& XS() const{ return *xs; }

  bool IsFaceLocal(int f) const {return face_local_flags[f];}

  int NumNodes() const {return num_nodes;}

  double Volume() const {return volume;}

  void ZeroOutflow(     ) {outflow.assign(outflow.size(),0.0);}
  void ZeroOutflow(int g) {if (g<outflow.size()) outflow[g]=0.0;}
  void AddOutflow(int g, double intS_mu_psi)
  {
    if (g<outflow.size()) outflow[g] += intS_mu_psi;
  }
  double GetOutflow(int g) const
  {
    if (g<outflow.size()) return outflow[g];
    else return 0.0;
  }

  void ReassingXS(const chi_physics::MultiGroupXSBase& xs_mapped)
  {
    xs = &xs_mapped;
  }
};


struct UnitCellMatrices
{
  MatDbl  K_matrix;
  MatVec3 G_matrix;
  MatDbl  M_matrix;
  VecDbl  Vi_vectors;

  std::vector<MatDbl>  face_M_matrices;
  std::vector<MatVec3> face_G_matrices;
  std::vector<VecDbl>  face_Si_vectors;
};

enum class AGSSchemeEntryType
{
  GROUPSET_ID = 1,
  SCHEME      = 2
};

class AGSSchemeEntry
{
private:
  const AGSSchemeEntryType type_;
  const int groupset_id_ = 0;
  const std::string scheme_name_;
  std::vector<AGSSchemeEntry> scheme_entries_;
public:
  explicit
  AGSSchemeEntry(int groupset_id) :
    type_(AGSSchemeEntryType::GROUPSET_ID),
    groupset_id_(groupset_id)
  {}

  explicit
  AGSSchemeEntry(const std::string& scheme) :
    type_(AGSSchemeEntryType::SCHEME),
    scheme_name_(scheme)
  {}

  AGSSchemeEntryType Type() const {return type_;}
  int GroupsetID() const {return groupset_id_;}

  std::vector<AGSSchemeEntry>& SchemeEntries() {return scheme_entries_;}
};

}//namespace lbs

#endif
