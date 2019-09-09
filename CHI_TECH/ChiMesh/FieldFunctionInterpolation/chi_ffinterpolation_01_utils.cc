#include "chi_ffinterpolation.h"
#include "../MeshHandler/chi_meshhandler.h"
#include "../VolumeMesher/chi_volumemesher.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/** This functions checks the intersection of a plane with a tetrahedron.
 * The equation of a plane is
 *      nx(x-x0) + ny(y-y0) + nz(z-z0) = 0
 * Where the plane normal is (nx,ny,nz) and the plane point is (x0,y0,z0).
 * If we form a dot product between the normal and a vector
 * (x-x0,y-y0,z-z0) then sign of the result gives the sense to the surface.
 * Therefore, if we encounter differing senses then the plane is indeed
 * intersecting.*/
bool chi_mesh::FieldFunctionInterpolation::
CheckPlaneTetIntersect(chi_mesh::Normal plane_normal,
                       chi_mesh::Vector plane_point,
                       std::vector<chi_mesh::Vector *>* tet_points)
{
  bool current_sense = false;

  size_t num_points = tet_points->size();
  for (int i=0; i<num_points; i++)
  {
    chi_mesh::Vector v = (*(*tet_points)[i]) - plane_point;
    double dotp = plane_normal.Dot(v);

    bool new_sense = (dotp >= 0.0);

    if (i==0)
      current_sense = new_sense;
    else if (new_sense != current_sense)
      return true;
  }
  return false;
}


//###################################################################
bool chi_mesh::FieldFunctionInterpolation::
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

//###################################################################
/**Equation of a plane:    nx*(x-x0) + ny*(y-y0) + nz*(z-z0) = 0
                   or      ax + by + cz = D    <a,b,c>=<nx,ny,nz>
                                                  D=<nx*x0 + ny*y0 + nz*z0
Equation of a line:     <x,y,z> = <x0,y0,z0> + d<x1,y1,z1>
                                = v0 + d*v1
                        E=n-dot

*/
bool chi_mesh::FieldFunctionInterpolation::
     CheckLineTriangleIntersect(std::vector<chi_mesh::Vector>& triangle_points,
                                chi_mesh::Vector line_point_i,
                                chi_mesh::Vector line_point_f)
{
  //======================================== First find plane intersection
  //Compute normal
  chi_mesh::Vector p01 = triangle_points[1]-triangle_points[0];
  chi_mesh::Vector p12 = triangle_points[2]-triangle_points[1];
  chi_mesh::Vector p20 = triangle_points[0]-triangle_points[2];
  chi_mesh::Vector pintersection;
  chi_mesh::Vector n   = p01.Cross(p12); n=n/n.Norm();

  std::pair<double,double> weights;
  if (CheckPlaneLineIntersect(n,triangle_points[0],
                              line_point_i,
                              line_point_f,
                              pintersection,weights))
  {
    //======================================== Now determine if it is inside
    //                                         the triangle
    chi_mesh::Vector p0_pint = pintersection - triangle_points[0];
    chi_mesh::Vector p1_pint = pintersection - triangle_points[1];
    chi_mesh::Vector p2_pint = pintersection - triangle_points[2];

    chi_mesh::Vector cross0 = p01.Cross(p0_pint); cross0=cross0/cross0.Norm();
    chi_mesh::Vector cross1 = p12.Cross(p1_pint); cross1=cross1/cross1.Norm();
    chi_mesh::Vector cross2 = p20.Cross(p2_pint); cross2=cross2/cross2.Norm();

    if (cross0.Dot(n)<(0.0)) return false;
    if (cross1.Dot(n)<(0.0)) return false;
    if (cross2.Dot(n)<(0.0)) return false;

    return true;
  }

  return false;
}



void chi_mesh::FieldFunctionInterpolation::
CreateCFEMMapping(int num_grps, int num_moms, int g, int m,
                  Vec& x, Vec& x_cell,
                  std::vector<int> &cfem_nodes,
                  std::vector<int> *mapping)
{
  chi_mesh::MeshHandler* cur_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher* mesher = cur_handler->volume_mesher;

  size_t num_nodes_to_map = cfem_nodes.size();
  std::vector<int> mapped_nodes;
  for (int n=0; n< num_nodes_to_map; n++)
  {
    int ir = mesher->MapNode(
      cfem_nodes[n])*num_grps + g;

    mapped_nodes.push_back(ir);
    mapping->push_back(n);
  }


  VecCreateSeq(PETSC_COMM_SELF,mapped_nodes.size()+1,&x_cell);
  VecSet(x_cell,0.0);


  IS global_set;
  IS local_set;
  ISCreateGeneral(PETSC_COMM_WORLD, num_nodes_to_map, mapped_nodes.data(),
                  PETSC_COPY_VALUES,&global_set);
  ISCreateGeneral(PETSC_COMM_WORLD, num_nodes_to_map, mapping->data(),
                  PETSC_COPY_VALUES,&local_set);
  VecScatter scat;
  VecScatterCreate(x,global_set,x_cell,local_set,&scat);
  VecScatterBegin(scat,x,x_cell,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(scat,x,x_cell,INSERT_VALUES,SCATTER_FORWARD);

  ISDestroy(&global_set);
  ISDestroy(&local_set);
}


//###################################################################
/**Computes interpolated field function values.*/
void chi_mesh::FieldFunctionInterpolation::
CreatePWLDMapping(int num_grps, int num_moms, int g, int m,
                  std::vector<int> &pwld_nodes,
                  std::vector<int> &pwld_cells,
                  std::vector<int> &local_cell_dof_array_address,
                  std::vector<int> *mapping)
{
  size_t num_nodes_to_map = pwld_nodes.size();
  for (int n=0; n< num_nodes_to_map; n++)
  {
    int dof = pwld_nodes[n];
    int cell_g_index = pwld_cells[n];

    auto cell = grid_view->cells[cell_g_index];

    int c = cell->cell_local_id;
    int dof_map_start = local_cell_dof_array_address[c];

    int address = dof_map_start + dof*num_grps*num_moms + num_grps*m + g;

    mapping->push_back(address);
  }

}

//###################################################################
/**Computes interpolated field function values.*/
void chi_mesh::FieldFunctionInterpolation::
  CreateFVMapping(int num_grps, int num_moms, int g, int m,
                  std::vector<int> &cells, std::vector<int> *mapping)
{
  for (int c=0; c<cells.size(); c++)
  {
    int cell_glob_index = cells[c];
    int address = cell_glob_index*num_grps*num_moms + num_grps*m + g;

    mapping->push_back(address);
  }
}