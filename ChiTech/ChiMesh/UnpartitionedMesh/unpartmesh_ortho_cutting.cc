/**The routines in this file are intended to handle orthogonal cuts of
 * a mesh for purposes of KBA partitioning. To make it simple the routines
 * only operate on unpartitioned meshes.
 *
 * The designs specifications here is to support both arbitrary orthogonal
 * cuts as well load-balanced by dimension cuts. In the case of the latter
 * we first divide the domain with a number of x-cuts. In-between these x-cuts
 * we have a sub-domain that can be subdivided in the y-direction independent
 * from other subdomains. In both the arbitrary and by-dimension case we assume
 * extruded-type geometry in the z-direction.
 *
 * The simplest mechanism with which we can cut a mesh is by using a plane.
 * This will cut across the entire domain and is useful for the x- and z-cuts.
 * A plane simply requires a reference point and a normal to complete its
 * definition. A cell will be
 * cut by a plane if its vertices appear on both sides of the plane.
 * The y-cuts, however, may be done by-dimension and therefore we need to
 * consider not just a simple plane but rather a strip. To complete the
 * definition of a strip we need two reference points and a normal. A cell is
 * intersected by a strip if any edge (forming a 3D-line) intersects the strip.
 * Since the logic of a line intersection a strip is used in the raytracing of
 * 2D cells, this functionality already exists in Chi-Tech.
 *
 * Cutting the mesh is normally just a consideration in 2D and 3D so we will
 * maintain our focus on these cases for now.
 *
 * */
#include "chi_unpartitioned_mesh.h"

//###################################################################
/**

*/