#include "pwl_polygon.h"

//###################################################################
/**Returns the value of the shape function given cartesian
 * coordinates.*/
double PolygonFEView::Shape_xy(int i, chi_mesh::Vector xyz)
{
  for (int s=0; s<num_of_subtris; s++)
  {
    chi_mesh::Vector p0 = *grid->nodes[sides[s]->v_index[0]];
    chi_mesh::Vector xyz_ref = xyz - p0;

    chi_mesh::Vector xi_eta_zeta   = sides[s]->Jinv*xyz_ref;

    double xi  = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    //Determine if inside tet
    if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and
        ((xi + eta)<=(1.0+1.0e-12)))
    {
      int index = node_to_side_map[i][s];
      double value = 0.0;

      if (index==0)
      {
        value = 1.0 - xi - eta;
      }
      if (index==1)
      {
        value = xi;
      }

      value += beta*eta;

      return value;
    }
  }

  return 0.0;
}

//###################################################################
/**Returns the value of the gradient of the shape function.*/
chi_mesh::Vector PolygonFEView::GradShape_xy(int i, chi_mesh::Vector xyz)
{
  chi_mesh::Vector grad_r;
  chi_mesh::Vector grad;

  for (int e=0; e<num_of_subtris; e++)
  {
    chi_mesh::Vector p0 = *grid->nodes[sides[e]->v_index[0]];
    chi_mesh::Vector xyz_ref = xyz - p0;

    chi_mesh::Vector xi_eta_zeta = sides[e]->Jinv*xyz_ref;

    double xi  = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;
    double zeta= xi_eta_zeta.z;

    if ((xi>=-1.0e-12) and (eta>=-1.0e-12)  and
        ((xi + eta)<=(1.0+1.0e-12)))
    {
      int index = node_to_side_map[i][e];

      if (index == 0)
      {
        grad_r.x += -1.0;
        grad_r.y += -1.0;
      }
      if (index == 1)
      {
        grad_r.x += 1.0;
        grad_r.y += 0.0;
      }

      grad_r.y += beta*1.0;

      grad = sides[e]->JTinv*grad_r;

      return grad;
    }
  }

  return grad;
}