#include "pwl_polygon.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Returns the value of the shape function given cartesian
 * coordinates.*/
double PolygonMappingFE_PWL::ShapeValue(const int i, const chi_mesh::Vector3& xyz)
{
  for (int s=0; s<num_of_subtris; s++)
  {
    const auto& p0 = grid->vertices[sides[s].v_index[0]];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 xi_eta_zeta   = sides[s].Jinv * xyz_ref;

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
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void PolygonMappingFE_PWL::ShapeValues(const chi_mesh::Vector3 &xyz,
                                       std::vector<double> &shape_values)
{
  shape_values.resize(num_nodes, 0.0);
  for (int s=0; s<num_of_subtris; s++)
  {
    const auto& p0 = grid->vertices[sides[s].v_index[0]];
    chi_mesh::Vector3 xi_eta_zeta   = sides[s].Jinv * (xyz - p0);

    double xi  = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    //Determine if inside tet
    if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and
        ((xi + eta)<=(1.0+1.0e-12)))
    {
      for (int i=0; i < num_nodes; i++)
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

        shape_values[i] = value;
      }
      return;
    }//if in triangle
  }//for side
}

//###################################################################
/**Returns the evaluation of grad-shape function i at the supplied point.*/
chi_mesh::Vector3 PolygonMappingFE_PWL::GradShapeValue(const int i,
                                                       const chi_mesh::Vector3& xyz)
{
  chi_mesh::Vector3 grad_r;
  chi_mesh::Vector3 grad;

  for (int e=0; e<num_of_subtris; e++)
  {
    const auto& p0 = grid->vertices[sides[e].v_index[0]];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 xi_eta_zeta = sides[e].Jinv * xyz_ref;

    double xi  = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

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

      grad = sides[e].JTinv*grad_r;

      return grad;
    }
  }

  return grad;
}

//###################################################################
/**Populates gradshape_values with the value of each shape function's
 * gradient evaluated at the supplied point.*/
void PolygonMappingFE_PWL::GradShapeValues(
  const chi_mesh::Vector3 &xyz,
  std::vector<chi_mesh::Vector3> &gradshape_values)
{
  gradshape_values.clear();
  for (int i=0; i < num_nodes; ++i)
    gradshape_values.emplace_back(GradShapeValue(i,xyz));
}