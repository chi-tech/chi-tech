#ifndef LBS_STRUCTS_H
#define LBS_STRUCTS_H

#include "math/chi_math.h"
#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"
#include "physics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include <functional>
#include <map>

namespace lbs
{

typedef std::vector<size_t> DirIDs; ///< Direction-IDs
typedef std::vector<DirIDs> UniqueSOGroupings;
typedef std::map<size_t, size_t> DirIDToSOMap;

typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> MatDbl;
typedef std::vector<chi_mesh::Vector3> VecVec3;
typedef std::vector<VecVec3> MatVec3;

typedef std::shared_ptr<chi_physics::MultiGroupXS> XSPtr;
typedef std::shared_ptr<chi_physics::IsotropicMultiGrpSource> IsotropicSrcPtr;

enum class SolverType
{
  DISCRETE_ORDINATES = 1,
  DIFFUSION_DFEM = 2,
  DIFFUSION_CFEM = 3,
};

enum class GeometryType
{
  NO_GEOMETRY_SET = 0,
  ONED_SLAB = 1,
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL = 3,
  TWOD_CARTESIAN = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6
};

inline chi_math::CoordinateSystemType
MapGeometryTypeToCoordSys(const GeometryType gtype)
{
  using namespace chi_math;
  switch (gtype)
  {
    case GeometryType::ONED_SLAB:
    case GeometryType::TWOD_CARTESIAN:
    case GeometryType::THREED_CARTESIAN:
      return CoordinateSystemType::CARTESIAN;
    case GeometryType::ONED_SPHERICAL:
      return CoordinateSystemType::SPHERICAL;
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
      return CoordinateSystemType::CYLINDRICAL;
    default:
      return CoordinateSystemType::CARTESIAN;
  }
}

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
  REFLECTING = 3,
  INCIDENT_ANISTROPIC_HETEROGENEOUS = 4
};

struct BoundaryPreference
{
  BoundaryType type;
  std::vector<double> isotropic_mg_source;
  std::string source_function;
};

enum SourceFlags : int
{
  NO_FLAGS_SET = 0,
  APPLY_FIXED_SOURCES = (1 << 0),
  APPLY_WGS_SCATTER_SOURCES = (1 << 1),
  APPLY_AGS_SCATTER_SOURCES = (1 << 2),
  APPLY_WGS_FISSION_SOURCES = (1 << 3),
  APPLY_AGS_FISSION_SOURCES = (1 << 4),
  SUPPRESS_WG_SCATTER = (1 << 5),
  ZERO_INCOMING_DELAYED_PSI = (1 << 6)
};

inline SourceFlags operator|(const SourceFlags f1, const SourceFlags f2)
{
  return static_cast<SourceFlags>(static_cast<int>(f1) | static_cast<int>(f2));
}

enum class PhiSTLOption
{
  PHI_OLD = 1,
  PHI_NEW = 2
};

class LBSGroupset;
typedef std::function<void(LBSGroupset& groupset,
                           std::vector<double>& destination_q,
                           const std::vector<double>& phi,
                           SourceFlags source_flags)>
  SetSourceFunction;

class AGSSchemeEntry;

/**Struct for storing LBS options.*/
struct Options
{
  typedef chi_math::SpatialDiscretizationType SDMType;

  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  SDMType sd_type = SDMType::PIECEWISE_LINEAR_DISCONTINUOUS;
  unsigned int scattering_order = 1;
  int sweep_eager_limit = 32000; // see chiLBSSetProperty documentation

  bool read_restart_data = false;
  std::string read_restart_folder_name = std::string("YRestart");
  std::string read_restart_file_base = std::string("restart");

  bool write_restart_data = false;
  std::string write_restart_folder_name = std::string("YRestart");
  std::string write_restart_file_base = std::string("restart");
  double write_restart_interval = 30.0;

  bool use_precursors = false;
  bool use_src_moments = false;

  bool save_angular_flux = false;

  bool verbose_inner_iterations = true;
  bool verbose_ags_iterations = false;
  bool verbose_outer_iterations = true;

  bool power_field_function_on = false;
  double power_default_kappa = 3.20435e-11; // 200MeV to Joule
  double power_normalization = -1.0;

  std::string field_function_prefix_option = "prefix";
  std::string field_function_prefix; // Default is empty

  Options() = default;

  std::vector<AGSSchemeEntry> ags_scheme;
};

/**Transport view of a cell*/
class CellLBSView
{
private:
  size_t phi_address_;
  int num_nodes_;
  int num_groups_;
  int num_grps_moms_;
  const chi_physics::MultiGroupXS* xs_;
  double volume_;
  const std::vector<bool> face_local_flags_;
  const std::vector<int> face_locality_;
  const std::vector<const chi_mesh::Cell*> neighbor_cell_ptrs_;
  std::vector<double> outflow_;

public:
  CellLBSView(size_t phi_address,
              int num_nodes,
              int num_groups,
              int num_moments,
              const chi_physics::MultiGroupXS& xs_mapping,
              double volume,
              const std::vector<bool>& face_local_flags,
              const std::vector<int>& face_locality,
              const std::vector<const chi_mesh::Cell*>& neighbor_cell_ptrs,
              bool cell_on_boundary)
    : phi_address_(phi_address),
      num_nodes_(num_nodes),
      num_groups_(num_groups),
      num_grps_moms_(num_groups * num_moments),
      xs_(&xs_mapping),
      volume_(volume),
      face_local_flags_(face_local_flags),
      face_locality_(face_locality),
      neighbor_cell_ptrs_(neighbor_cell_ptrs)
  {
    if (cell_on_boundary) outflow_.resize(num_groups_, 0.0);
  }

  size_t MapDOF(int node, int moment, int grp) const
  {
    return phi_address_ + node * num_grps_moms_ + num_groups_ * moment + grp;
  }

  const chi_physics::MultiGroupXS& XS() const { return *xs_; }

  bool IsFaceLocal(int f) const { return face_local_flags_[f]; }
  int FaceLocality(int f) const { return face_locality_[f]; }
  const chi_mesh::Cell* FaceNeighbor(int f) const
  {
    return neighbor_cell_ptrs_[f];
  }

  int NumNodes() const { return num_nodes_; }

  double Volume() const { return volume_; }

  void ZeroOutflow() { outflow_.assign(outflow_.size(), 0.0); }
  void ZeroOutflow(int g)
  {
    if (g < outflow_.size()) outflow_[g] = 0.0;
  }
  void AddOutflow(int g, double intS_mu_psi)
  {
    if (g < outflow_.size()) outflow_[g] += intS_mu_psi;
  }

  double GetOutflow(int g) const
  {
    if (g < outflow_.size()) return outflow_[g];
    else
      return 0.0;
  }

  void ReassingXS(const chi_physics::MultiGroupXS& xs_mapped)
  {
    xs_ = &xs_mapped;
  }
};

struct UnitCellMatrices
{
  MatDbl K_matrix;
  MatVec3 G_matrix;
  MatDbl M_matrix;
  VecDbl Vi_vectors;

  std::vector<MatDbl> face_M_matrices;
  std::vector<MatVec3> face_G_matrices;
  std::vector<VecDbl> face_Si_vectors;
};

enum class AGSSchemeEntryType
{
  GROUPSET_ID = 1,
  SCHEME = 2
};

class AGSSchemeEntry
{
private:
  const AGSSchemeEntryType type_;
  const int groupset_id_ = 0;
  const std::string scheme_name_;
  std::vector<AGSSchemeEntry> scheme_entries_;

public:
  explicit AGSSchemeEntry(int groupset_id)
    : type_(AGSSchemeEntryType::GROUPSET_ID), groupset_id_(groupset_id)
  {
  }

  explicit AGSSchemeEntry(const std::string& scheme)
    : type_(AGSSchemeEntryType::SCHEME), scheme_name_(scheme)
  {
  }

  AGSSchemeEntryType Type() const { return type_; }
  int GroupsetID() const { return groupset_id_; }

  std::vector<AGSSchemeEntry>& SchemeEntries() { return scheme_entries_; }
};

} // namespace lbs

#endif
