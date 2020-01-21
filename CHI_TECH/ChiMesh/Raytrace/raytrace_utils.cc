#include "raytracing.h"
#include <ChiMesh/Cell/cell.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include <algorithm>

//###################################################################
/**Computes the intersection of a line with a plane.
 *
 * The first step of this algorithm is to compute v0 and v1. These
 * are vectors from the plane's reference point to each of the line-points,
 * respectively.
 * We then take the dot-products of these
 * vectors with the plane normal.
 * We then say that the vectors have a positive sense if the dot-product
 * is positive and a negative sense if the dot-product is negative.
 * If the senses are not equal then the line intersects the plane.
 *
 * Since the face normal is a normalized vector the dot-product of v0 or v1
 * will give the projection of the relevant vector along the normal to the
 * plane. We can use this projection to compute a weight associated with
 * each vector. This also then allows us to compute the intersection point.

 \param plane_normal The normal associated with the plane
 \param plane_point The reference point for the plane
 \param line_point_0 The line's initial point
 \param line_point_1 The line's destination point
 \param intersection_point The point to be populated with the intersection
                           point
 \param weights The weights associated with this intersection

 \return Returns true if the line intersects the plane and false otherwise.

 \author Jan*/
bool chi_mesh::
CheckPlaneLineIntersect(const chi_mesh::Normal& plane_normal,
                        const chi_mesh::Vector& plane_point,
                        const chi_mesh::Vector& line_point_0,
                        const chi_mesh::Vector& line_point_1,
                        chi_mesh::Vector& intersection_point,
                        std::pair<double,double>& weights)
{
  chi_mesh::Vector v0 = line_point_0 - plane_point;
  chi_mesh::Vector v1 = line_point_1 - plane_point;

  double dotp_0 = plane_normal.Dot(v0);
  double dotp_1 = plane_normal.Dot(v1);

  bool sense_0 = (dotp_0 >= 0.0);
  bool sense_1 = (dotp_1 >= 0.0);

  if (sense_0 != sense_1)
  {
    double dotp_total = std::fabs(dotp_0) + std::fabs(dotp_1);
    weights.first = (std::fabs(dotp_0)/dotp_total);
    weights.second = 1.0 - weights.first;
    intersection_point =
      line_point_0*weights.second +
      line_point_1*weights.first;

    return true;
  }

  return false;
}

//###################################################################
/**Given a strip defined by two points (v0,v1) and a normal, n,
 * (meaning infinite in the direction defined by (v1-v0).cross(n),
 * this function determines if a line, defined from p0 to p1,
 * intersects it. If it does then `true` is returned and
 * `intersection_point` contains the point of intersection. If it does
 * not then `false` is returned and `intersection_point` remains
 * unchanged.
 *
 * */
bool chi_mesh::CheckLineIntersectStrip(
  const chi_mesh::Vector& strip_point0,
  const chi_mesh::Vector& strip_point1,
  const chi_mesh::Vector& strip_normal,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  chi_mesh::Vector& intersection_point)
{
  chi_mesh::Vector plane_intersection_point;
  std::pair<double,double> weights;

  bool intersects_plane = chi_mesh::CheckPlaneLineIntersect(
    strip_normal, strip_point0,
    line_point0, line_point1,
    plane_intersection_point, weights);

  if (!intersects_plane) return false;

  chi_mesh::Vector edge_vec = strip_point1 - strip_point0;
  chi_mesh::Vector ints_vec1 = plane_intersection_point - strip_point0;
  chi_mesh::Vector ints_vec2 = plane_intersection_point - strip_point1;

  bool sense1 = edge_vec.Dot(ints_vec1)>=0.0;
  bool sense2 = edge_vec.Dot(ints_vec2)>=0.0;

  if (sense1 != sense2)
  {
    intersection_point = plane_intersection_point;

    return true;
  }

  return false;
}

//###################################################################
/**Given a triangle defined by three points, computes whether a line
 * intersects this triangle.
 *
 * */
bool chi_mesh::CheckLineIntersectTriangle(
  const chi_mesh::Vector& tri_point0,
  const chi_mesh::Vector& tri_point1,
  const chi_mesh::Vector& tri_point2,
  const chi_mesh::Normal& tri_normal,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  chi_mesh::Vector& intersection_point)
{
  chi_mesh::Vector plane_pi; //plane point of intersection
  std::pair<double,double> weights;

  bool intersects_plane =
    chi_mesh::CheckPlaneLineIntersect(
      tri_normal, tri_point0,
      line_point0, line_point1,
      plane_pi, weights);

  if (!intersects_plane) return false;
  else
  {
    bool in_triangle =
      chi_mesh::CheckPointInTriangle(
        tri_point0, tri_point1, tri_point2, tri_normal, plane_pi);

    if (in_triangle)
    {
      intersection_point = plane_pi;
      return true;
    }
  }

  return false;
}

bool
chi_mesh::CheckLineIntersectTriangle2(
  const chi_mesh::Vector& tri_point0,
  const chi_mesh::Vector& tri_point1,
  const chi_mesh::Vector& tri_point2,
  const chi_mesh::Vector& ray_posi,
  const chi_mesh::Vector& ray_dir,
  chi_mesh::Vector& intersection_point)
{
  double epsilon = 1.0e-12;
  chi_mesh::Vector edge1 = tri_point1 - tri_point0;
  chi_mesh::Vector edge2 = tri_point2 - tri_point0;

  // Compute characteristic vector for incident angle
  // This vector becomes perpendicular to the plane
  // when the ray is parallel to triangle
  chi_mesh::Vector h = ray_dir.Cross(edge2);

  // If h is indeed perpendicular to the plane,
  // the dot product of the other leg with this h will be close
  // to zero.
  double a = edge1.Dot(h);
  if (std::fabs(a) < epsilon)
    return false;

  chi_mesh::Vector s = ray_posi - tri_point0;

  double f = 1.0/a;

  // If, s projected onto h, is greater than,
  // v01 projected onto h, or negative, there is now way
  // the ray can intersect
  double u = f*(s.Dot(h));
  if (u<0.0 or u>1.0)
    return false;

  chi_mesh::Vector q = s.Cross(edge1);

  // If, q projected onto omega, is greater than,
  // v01 projected onto h, or negative, there is now way
  // the ray can intersect
  double v = f*ray_dir.Dot(q);
  if (v<0.0 or (u+v)>1.0)
    return false;

  double t = f*edge2.Dot(q);
  if (t > epsilon and t<(1.0/epsilon))
  {
    intersection_point = ray_posi + ray_dir*t;
    return true;
  } else
  {
    return false;
  }

}

//###################################################################
/** Check whether a point lies in a triangle.*/
bool
chi_mesh::CheckPointInTriangle(
  const chi_mesh::Vector& v0,
  const chi_mesh::Vector& v1,
  const chi_mesh::Vector& v2,
  const chi_mesh::Normal& n,
  const chi_mesh::Vector& point)
{
  auto v01 = v1 - v0;
  auto v12 = v2 - v1;
  auto v20 = v0 - v2;

  auto v0p = point - v0;
  auto v1p = point - v1;
  auto v2p = point - v2;

  auto vc0 = v01.Cross(v0p);
  auto vc1 = v12.Cross(v1p);
  auto vc2 = v20.Cross(v2p);

  bool dp0 = (vc0.Dot(n) >= 0.0);
  bool dp1 = (vc1.Dot(n) >= 0.0);
  bool dp2 = (vc2.Dot(n) >= 0.0);

  if (dp0 and dp1 and dp2)
    return true;
  else
    return false;

}

//###################################################################
/** Populates segment lengths along a ray.*/
void chi_mesh::PopulateRaySegmentLengths(
  const chi_mesh::MeshContinuum* grid,
  Cell* cell,
  std::vector<double> &segment_lengths,
  const chi_mesh::Vector& line_point0,
  const chi_mesh::Vector& line_point1,
  const chi_mesh::Vector& omega)
{
  double track_length = (line_point1-line_point0).Norm();
  if (segment_lengths.empty())
    segment_lengths.push_back(track_length);

  //======================================== Determine intersection points
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  // Since there are no segments within a slab we will only have
  // a single segment length. It is already pushed

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  // A polygon can be decomposed into "sides" by means of its
  // edges. Each side comprises a triangle formed by: the two
  // vertices of the associated edge v0 and v1, and the cell
  // centroid vc.
  // Since the triangles all share an edge we only determine
  // segment lengths from the strip defined by v0 to vc.
  if (cell->Type() == chi_mesh::CellType::POLYGON)
  {
    auto poly_cell = (chi_mesh::CellPolygon*)cell;

    int f=-1;
    for (auto& face : cell->faces) //edges
    {
      f++;
      chi_mesh::Vertex& v0 = *grid->nodes[face.vertex_ids[0]];
      chi_mesh::Vertex& vc = cell->centroid;

      auto& n0 = poly_cell->GetSegmentNormals(grid)[f];

      chi_mesh::Vertex intersection_point;
      bool intersects = chi_mesh::CheckLineIntersectStrip(
        v0, vc, n0,
        line_point0, line_point1,
        intersection_point);

      if (intersects)
      {
        double d = (intersection_point - line_point0).Norm();
        segment_lengths.push_back(d);
      }

    }//for face
  }
  else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    auto& vcc = cell->centroid;

    int f=-1;
    for (auto& face : cell->faces)
    {
      f++;
      auto& vfc  = face.centroid;

      int s=-1;
      for (auto vi : face.vertex_ids)
      {
        s++;
        auto& vert = *grid->nodes[vi];

        chi_mesh::Vertex intersection_point;

        bool intersects = chi_mesh::CheckLineIntersectTriangle2(
          vert,vfc,vcc,line_point0,omega,intersection_point);

        if (intersects)
        {
          double d = (intersection_point - line_point0).Norm();
          if (d < track_length)
            segment_lengths.push_back(d);
        }
      }//for edge
    }//for face
  }

  //======================================== Sort distances
  std::sort(segment_lengths.begin(),
            segment_lengths.end());

  //======================================== Populate segment lengths
  //if there are N segments intersected then there will always be
  //N+1 distances.
  size_t num_distances = segment_lengths.size();
  double last_distance = segment_lengths[0];
  for (size_t di=1; di < num_distances; di++)
  {
    double new_seg_length =
      segment_lengths[di] - last_distance;
    last_distance = segment_lengths[di];
    segment_lengths[di] = new_seg_length;
  }

}