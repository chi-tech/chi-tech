#ifndef _chi_mesh_raytracing_h
#define _chi_mesh_raytracing_h

#include "../chi_mesh.h"

namespace chi_mesh
{

struct RayDestinationInfo
{
  int destination_face_neighbor; //Neighbor cell global-id
  int destination_face_index;    //Neighbor cell face index (deprecated)

  RayDestinationInfo() :
  destination_face_neighbor(-1),
  destination_face_index(-1)
  {}
};

//=================================== Raytracing
RayDestinationInfo RayTrace(chi_mesh::MeshContinuum* grid,
                            Cell* cell,
                            const Vector3& pos_i,
                            const Vector3& omega_i,
                            double& d_to_surface,
                            Vector3& pos_f,
                            double epsilon_nudge=1.0e-8,
                            double backward_tolerance=1.0e-10,
                            double extension_distance=1.0e5,
                            int func_depth=0);

bool
CheckPlaneLineIntersect(const chi_mesh::Normal& plane_normal,
                        const chi_mesh::Vector3& plane_point,
                        const chi_mesh::Vector3& line_point_0,
                        const chi_mesh::Vector3& line_point_1,
                        chi_mesh::Vector3& intersection_point,
                        std::pair<double,double>& weights);

bool
CheckLineIntersectStrip(
  const chi_mesh::Vector3& strip_point0,
  const chi_mesh::Vector3& strip_point1,
  const chi_mesh::Vector3& strip_normal,
  const chi_mesh::Vector3& line_point0,
  const chi_mesh::Vector3& line_point1,
  chi_mesh::Vector3& intersection_point,
  double* distance_to_intersection = nullptr);

bool
CheckLineIntersectTriangle2(
  const chi_mesh::Vector3& tri_point0,
  const chi_mesh::Vector3& tri_point1,
  const chi_mesh::Vector3& tri_point2,
  const chi_mesh::Vector3& ray_posi,
  const chi_mesh::Vector3& ray_dir,
  chi_mesh::Vector3& intersection_point,
  double* distance_to_intersection = nullptr);

void PopulateRaySegmentLengths(
  const chi_mesh::MeshContinuum* grid,
  Cell* cell,
  std::vector<double> &segment_lengths,
  const chi_mesh::Vector3& line_point0,
  const chi_mesh::Vector3& line_point1,
  const chi_mesh::Vector3& omega);

bool
CheckPointInTriangle(const chi_mesh::Vector3& v0,
                     const chi_mesh::Vector3& v1,
                     const chi_mesh::Vector3& v2,
                     const chi_mesh::Normal& n,
                     const chi_mesh::Vector3& point);

}
#endif