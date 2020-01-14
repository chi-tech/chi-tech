#ifndef _chi_mesh_raytracing_h
#define _chi_mesh_raytracing_h

#include "../chi_mesh.h"

namespace chi_mesh
{

struct RayDestinationInfo
{
  int destination_face_neighbor;
  int destination_face_index;

  RayDestinationInfo() :
  destination_face_neighbor(-1),
  destination_face_index(-1)
  {}
};

//=================================== Raytracing
RayDestinationInfo RayTrace(chi_mesh::MeshContinuum* grid,
                            Cell* cell,
                            const Vector& pos_i,
                            const Vector& omega_i,
                            double& d_to_surface,
                            Vector& pos_f,
                            int func_depth=0);

bool
CheckPlaneLineIntersect(const chi_mesh::Normal& plane_normal,
                        const chi_mesh::Vector& plane_point,
                        const chi_mesh::Vector& line_point_0,
                        const chi_mesh::Vector& line_point_1,
                        chi_mesh::Vector& intersection_point,
                        std::pair<double,double>& weights);

bool
CheckLineIntersectStrip(
  const chi_mesh::Vector& strip_point0,
  const chi_mesh::Vector& strip_point1,
  const chi_mesh::Vector& strip_normal,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  chi_mesh::Vector& intersection_point);

bool
CheckLineIntersectTriangle(
  const chi_mesh::Vector& tri_point0,
  const chi_mesh::Vector& tri_point1,
  const chi_mesh::Vector& tri_point2,
  const chi_mesh::Normal& tri_normal,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  chi_mesh::Vector& intersection_point);

bool
CheckLineIntersectTriangle2(
  const chi_mesh::Vector& tri_point0,
  const chi_mesh::Vector& tri_point1,
  const chi_mesh::Vector& tri_point2,
  const chi_mesh::Vector& ray_posi,
  const chi_mesh::Vector& ray_dir,
  chi_mesh::Vector& intersection_point);

void PopulateRaySegmentLengths(
  const chi_mesh::MeshContinuum* grid,
  Cell* cell,
  std::vector<double> &segment_lengths,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  const chi_mesh::Vector& omega);

bool
CheckPointInTriangle(const chi_mesh::Vector& v0,
                     const chi_mesh::Vector& v1,
                     const chi_mesh::Vector& v2,
                     const chi_mesh::Normal& n,
                     const chi_mesh::Vector& point);

}
#endif