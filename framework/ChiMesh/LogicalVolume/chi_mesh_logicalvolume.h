#ifndef CHI_MESH_LOGICALVOLUME_H
#define CHI_MESH_LOGICALVOLUME_H

#include "../chi_mesh.h"
#include <chi_log.h>
#include <array>

namespace chi_mesh
{
enum class LogicalVolumeType
{
  LVSPHERE = 1,
  LVSPHERE_ORIGIN = 2,
  LVRPP = 3,
  LVRCC = 4,
  LVSURFACE = 9,
  LVBOOLEAN = 10
};

// ###################################################################
/** Class for defining base logical volumes.*/
class LogicalVolume
{
private:
  const LogicalVolumeType type_;

protected:
  explicit LogicalVolume(LogicalVolumeType type) : type_(type) {}

public:
  LogicalVolumeType Type() const { return type_; }

  virtual bool Inside(const chi_mesh::Vector3& point) const { return false; }
};

} // namespace chi_mesh

// ###################################################################
/**Spherical logical volume.*/
class chi_mesh::SphereLogicalVolume : public LogicalVolume
{
public:
  double r_;
  double x0_, y0_, z0_;

  SphereLogicalVolume() : LogicalVolume(LogicalVolumeType::LVSPHERE_ORIGIN)
  {
    r_ = 1.0;
    x0_ = 0.0;
    y0_ = 0.0;
    z0_ = 0.0;
  }

  explicit SphereLogicalVolume(double in_radius)
    : LogicalVolume(LogicalVolumeType::LVSPHERE_ORIGIN)
  {
    r_ = in_radius;
    x0_ = 0.0;
    y0_ = 0.0;
    z0_ = 0.0;
  }

  SphereLogicalVolume(double in_radius, double in_x, double in_y, double in_z)
    : LogicalVolume(LogicalVolumeType::LVSPHERE)
  {
    r_ = in_radius;
    x0_ = in_x;
    y0_ = in_y;
    z0_ = in_z;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    double dx = point.x - x0_;
    double dy = point.y - y0_;
    double dz = point.z - z0_;

    double R2 = dx * dx + dy * dy + dz * dz;

    if (R2 <= (r_ * r_)) return true;
    else
      return false;
  }
};

// ###################################################################
/**Rectangular Parallel Piped (RPP) logical volume*/
class chi_mesh::RPPLogicalVolume : public LogicalVolume
{
public:
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;

  RPPLogicalVolume() : LogicalVolume(LogicalVolumeType::LVRPP)
  {
    xmin_ = 0.0;
    xmax_ = 1.0;
    ymin_ = 0.0;
    ymax_ = 1.0;
    zmin_ = 0.0;
    zmax_ = 1.0;
  }

  RPPLogicalVolume(
    double x0, double x1, double y0, double y1, double z0, double z1)
    : LogicalVolume(LogicalVolumeType::LVRPP)
  {
    xmin_ = x0;
    xmax_ = x1;
    ymin_ = y0;
    ymax_ = y1;
    zmin_ = z0;
    zmax_ = z1;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    if ((point.x <= xmax_) && (point.x >= xmin_) && (point.y <= ymax_) &&
        (point.y >= ymin_) && (point.z <= zmax_) && (point.z >= zmin_))
    {
      return true;
    }
    else
      return false;
  }
};

// ###################################################################
/**Right Circular Cylinder (RCC) logical volume.
 *
 * Determining whether a point is within an RCC is tricky.
 * */
class chi_mesh::RCCLogicalVolume : public LogicalVolume
{
public:
  double x0_, y0_, z0_;
  double vx_, vy_, vz_;
  double r_;

  RCCLogicalVolume() : LogicalVolume(LogicalVolumeType::LVRCC)
  {
    x0_ = 0.0;
    y0_ = 0.0;
    z0_ = 0.0;
    vx_ = 0.0;
    vy_ = 0.0;
    vz_ = 1.0;
    r_ = 1.0;
  }

  RCCLogicalVolume(double ix0,
                   double iy0,
                   double iz0,
                   double ivx,
                   double ivy,
                   double ivz,
                   double ir)
    : LogicalVolume(LogicalVolumeType::LVRCC)
  {
    x0_ = ix0;
    y0_ = iy0;
    z0_ = iz0;
    vx_ = ivx;
    vy_ = ivy;
    vz_ = ivz;
    r_ = ir;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    typedef chi_mesh::Vector3 Vec3;

    const auto& pr = point;                // reference point
    const Vec3 p0(x0_, y0_, z0_);          // cylinder root
    const Vec3 cyl_dir_vec(vx_, vy_, vz_); // cylinder direction vector
    const Vec3 k_hat(0.0, 0.0, 1.0);       // k_hat

    const Vec3 p0r = pr - p0;
    const Vec3 cyl_unit_dir = cyl_dir_vec.Normalized(); // aka cud
    const double cyl_length = cyl_dir_vec.Norm();

    //====================================== Check if point is within
    //                                       normal extents
    const double p0r_dot_cud = p0r.Dot(cyl_unit_dir);
    if (p0r_dot_cud < 0.0 or p0r_dot_cud > cyl_length) return false;

    //====================================== Building rotation matrix
    // This rotation matrix must be such that,
    // when a coordinate system is rotated with it,
    // the new normal vector points along the
    Vec3 binorm;
    Vec3 tangent;
    if (std::abs(cyl_dir_vec.Dot(k_hat) / cyl_dir_vec.Norm()) > (1.0 - 1.0e-12))
    {
      binorm = Vec3(0.0, 1.0, 0.0);
      tangent = Vec3(1.0, 0.0, 0.0);
    }
    else
    {
      binorm = k_hat.Cross(cyl_dir_vec);
      binorm = binorm / binorm.Norm();
      tangent = binorm.Cross(cyl_dir_vec);
      tangent = tangent / tangent.Norm();
    }

    //====================================== Project p0r onto the binorm and
    //                                       tangent
    const Vec3 p0r_projected(p0r.Dot(tangent), p0r.Dot(binorm), 0.0);

    //====================================== Determine if point is within
    //cylinder
    if (p0r_projected.NormSquare() <= r_ * r_) return true;
    else
      return false;
  }
};

// ###################################################################
/**SurfaceMesh volume*/
class chi_mesh::SurfaceMeshLogicalVolume : public LogicalVolume
{
private:
  typedef std::shared_ptr<chi_mesh::SurfaceMesh> SurfaceMeshPtr;
  std::array<double, 2> xbounds_;
  std::array<double, 2> ybounds_;
  std::array<double, 2> zbounds_;

public:
  const SurfaceMeshPtr surf_mesh = nullptr;

  explicit SurfaceMeshLogicalVolume(SurfaceMeshPtr in_surf_mesh);

  bool Inside(const chi_mesh::Vector3& point) const override;
};

// ###################################################################
/**Boolean volume*/
class chi_mesh::BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool, std::shared_ptr<LogicalVolume>>> parts;

  BooleanLogicalVolume() : LogicalVolume(LogicalVolumeType::LVBOOLEAN) {}

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    for (const auto& part : parts)
    {
      if (not(part.first && part.second->Inside(point))) return false;
    }
    return true;
  }
};

#endif // CHI_MESH_LOGICALVOLUME_H
