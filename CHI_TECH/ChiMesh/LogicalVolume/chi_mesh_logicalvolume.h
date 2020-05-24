#ifndef _chi_mesh_logicalvolume_h
#define _chi_mesh_logicalvolume_h

#include "../chi_mesh.h"
#include <chi_log.h>

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

  virtual bool Inside(chi_mesh::Vector3 point)
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

  SphereLogicalVolume(double in_radius)
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

  bool Inside(chi_mesh::Vector3 point)
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

  bool Inside(chi_mesh::Vector3 point)
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

  bool Inside(chi_mesh::Vector3 p1)
  {
    auto& chi_log = ChiLog::GetInstance();
    chi_mesh::Vector3 p0(x0, y0, z0);
    chi_mesh::Vector3 vd(vx, vy, vz);
    chi_mesh::Vector3 k(0.0, 0.0, 1.0);

    chi_mesh::Vector3 p01 = p1 - p0;

    //====================================== Building rotation matrix
    chi_mesh::Vector3 binorm;
    chi_mesh::Vector3 tangent;
    if (abs(vd.Dot(k)/vd.Norm())>(1.0-1.0e-12))
    {
      binorm = chi_mesh::Vector3(0.0, 1.0, 0.0);
      tangent = chi_mesh::Vector3(1.0, 0.0, 0.0);
    }
    else
    {
      binorm = k.Cross(vd);
      binorm = binorm/binorm.Norm();
      tangent = binorm.Cross(vd);
      tangent = tangent/tangent.Norm();
    }

    chi_mesh::Matrix3x3 R;

    R.SetColJVec(0,tangent);
    R.SetColJVec(1,binorm);
    R.SetColJVec(2,vd);

    //====================================== Rotate point to ref coords
    chi_mesh::Matrix3x3 Rinv = R.Inverse();
    chi_mesh::Vector3 p01T = Rinv * p01;

    chi_log.Log(LOG_0) << "Inverted p: \n" << p01T.PrintS();

    //====================================== Determine if point is within cylinder
    bool in_cylinder = false;
    double r2 = p01T.x*p01T.x + p01T.y*p01T.y;

    if (r2 <(this->r*this->r))
      in_cylinder = true;
    else
      return false;

    //====================================== Determine if point is within extents
    double dotP = p01T.Dot(vd);
    if ((dotP>=0.0) and (dotP<=1.0))
      return true;

    return false;
  }
};

//###################################################################
/**SurfaceMesh volume*/
class chi_mesh::SurfaceMeshLogicalVolume : public LogicalVolume
{
private:
  double xbounds[2];
  double ybounds[2];
  double zbounds[2];
public:
  chi_mesh::SurfaceMesh* surf_mesh;

  SurfaceMeshLogicalVolume(chi_mesh::SurfaceMesh* in_surf_mesh);

  bool Inside(chi_mesh::Vector3 point);
private:
  bool CheckPlaneLineIntersect(chi_mesh::Normal plane_normal,
                               chi_mesh::Vector3 plane_point,
                               chi_mesh::Vector3 line_point_0,
                               chi_mesh::Vector3 line_point_1,
                               chi_mesh::Vector3& intersection_point,
                               std::pair<double,double>& weights);
};


//###################################################################
/**Boolean volume*/
class chi_mesh::BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool,LogicalVolume*>> parts;

  bool Inside(chi_mesh::Vector3 point)
  {
    for (int p=0;p<parts.size();p++)
    {
      if (not (parts[p].first && parts[p].second->Inside(point)))
        return false;
    }
    return true;
  }
};


#endif