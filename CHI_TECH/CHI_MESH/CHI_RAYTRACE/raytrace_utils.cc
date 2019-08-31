#include "../chi_mesh.h"

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
CheckPlaneLineIntersect(chi_mesh::Normal plane_normal,
                        chi_mesh::Vector plane_point,
                        chi_mesh::Vector line_point_0,
                        chi_mesh::Vector line_point_1,
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