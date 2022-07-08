#ifndef _chi_mesh_logicalvolume_h
#define _chi_mesh_logicalvolume_h

#include "../chi_mesh.h"
#include <chi_log.h>
#include <array>

#define SPHERE        1
#define SPHERE_ORIGIN 2
#define RPP           3
#define RCC           4
#define SURFACE       9
#define BOOLEAN       10

//###################################################################
/** Class for defining base logical volumes.*/
class chi_mesh::LogicalVolume
{
public:
  int type_index;

  LogicalVolume()
  {
    type_index = 0;
  }

  int Type() const
  {
    return type_index;
  }

  virtual bool Inside(const chi_mesh::Vector3& point) const
  {
    return false;
  }
};

//###################################################################
/**Spherical logical volume.*/
class chi_mesh::SphereLogicalVolume : public LogicalVolume
{
public:
  double r;
  double x0,y0,z0;

  SphereLogicalVolume()
  {
    type_index = SPHERE_ORIGIN;
    r  = 1.0; x0 = 0.0; y0 = 0.0; z0 = 0.0;
  }

  explicit SphereLogicalVolume(double in_radius)
  {
    type_index = SPHERE_ORIGIN;
    r  = in_radius; x0 = 0.0; y0 = 0.0; z0 = 0.0;
  }

  SphereLogicalVolume(double in_radius,double in_x,double in_y, double in_z)
  {
    type_index = SPHERE;
    r  = in_radius;
    x0 = in_x;
    y0 = in_y;
    z0 = in_z;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    double dx = point.x - x0;
    double dy = point.y - y0;
    double dz = point.z - z0;

    double R2 = dx*dx + dy*dy + dz*dz;

    if (R2<= (r*r))
      return true;
    else
      return false;
  }
};

//###################################################################
/**Rectangular Parallel Piped (RPP) logical volume*/
class chi_mesh::RPPLogicalVolume : public LogicalVolume
{
public:
  double xmin,xmax;
  double ymin,ymax;
  double zmin,zmax;

  RPPLogicalVolume()
  {
    type_index = RPP;
    xmin = 0.0; xmax = 1.0;
    ymin = 0.0; ymax = 1.0;
    zmin = 0.0; zmax = 1.0;
  }

  RPPLogicalVolume(double x0, double x1,
                   double y0, double y1,
                   double z0, double z1)
  {
    type_index = RPP;
    xmin = x0; xmax = x1;
    ymin = y0; ymax = y1;
    zmin = z0; zmax = z1;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    if ((point.x <= xmax) && (point.x >= xmin) &&
        (point.y <= ymax) && (point.y >= ymin) &&
        (point.z <= zmax) && (point.z >= zmin))
    {
      return true;
    }
    else
      return false;
  }
};

//###################################################################
/**Right Circular Cylinder (RCC) logical volume.
 *
 * Determining whether a point is within an RCC is tricky.
 * */
class chi_mesh::RCCLogicalVolume : public LogicalVolume
{
public:
  double x0,y0,z0;
  double vx,vy,vz;
  double r;

  RCCLogicalVolume()
  {
    type_index = RCC;
    x0 = 0.0; y0 = 0.0; z0 = 0.0;
    vx = 0.0; vy = 0.0; vz = 1.0;
    r = 1.0;
  }

  RCCLogicalVolume(double ix0, double iy0, double iz0,
                   double ivx, double ivy, double ivz,
                   double ir)
  {
    type_index = RCC;
    x0 = ix0; y0 = iy0; z0 = iz0;
    vx = ivx; vy = ivy; vz = ivz;
    r = ir;
  }

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    typedef chi_mesh::Vector3 Vec3;
    auto& chi_log = ChiLog::GetInstance();

    const auto& pr = point;              //reference point
    const Vec3  p0(x0, y0, z0);          //cylinder root
    const Vec3  cyl_dir_vec(vx, vy, vz); //cylinder direction vector
    const Vec3  k_hat(0.0, 0.0, 1.0);    //k_hat

    const Vec3 p0r = pr - p0;
    const Vec3 cyl_unit_dir = cyl_dir_vec.Normalized(); //aka cud
    const double cyl_length = cyl_dir_vec.Norm();

    //====================================== Check if point is within
    //                                       normal extents
    const double p0r_dot_cud = p0r.Dot(cyl_unit_dir);
    if (p0r_dot_cud < 0.0 or p0r_dot_cud > cyl_length)
      return false;

    //====================================== Building rotation matrix
    //This rotation matrix must be such that,
    //when a coordinate system is rotated with it,
    //the new normal vector points along the
    Vec3 binorm;
    Vec3 tangent;
    if (std::abs(cyl_dir_vec.Dot(k_hat) / cyl_dir_vec.Norm()) > (1.0 - 1.0e-12))
    {
      binorm  = Vec3(0.0, 1.0, 0.0);
      tangent = Vec3(1.0, 0.0, 0.0);
    }
    else
    {
      binorm  = k_hat.Cross(cyl_dir_vec);
      binorm  = binorm/binorm.Norm();
      tangent = binorm.Cross(cyl_dir_vec);
      tangent = tangent/tangent.Norm();
    }

    //====================================== Project p0r onto the binorm and
    //                                       tangent
    const Vec3 p0r_projected(p0r.Dot(tangent), p0r.Dot(binorm), 0.0);

    //====================================== Determine if point is within cylinder
    if (p0r_projected.NormSquare() <= r*r)
      return true;
    else
      return false;
  }
};

//###################################################################
/**SurfaceMesh volume*/
class chi_mesh::SurfaceMeshLogicalVolume : public LogicalVolume
{
private:
  std::array<double,2> xbounds;
  std::array<double,2> ybounds;
  std::array<double,2> zbounds;
public:
  chi_mesh::SurfaceMesh* surf_mesh = nullptr;

  explicit
  SurfaceMeshLogicalVolume(chi_mesh::SurfaceMesh* in_surf_mesh);

  bool Inside(const chi_mesh::Vector3& point) const override;
};


//###################################################################
/**Boolean volume*/
class chi_mesh::BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool,LogicalVolume*>> parts;

  bool Inside(const chi_mesh::Vector3& point) const override
  {
    for (size_t p=0; p<parts.size();p++)
    {
      if (not (parts[p].first && parts[p].second->Inside(point)))
        return false;
    }
    return true;
  }
};


#endif
