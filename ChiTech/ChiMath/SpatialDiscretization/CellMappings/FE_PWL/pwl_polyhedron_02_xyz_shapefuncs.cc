#include "pwl_polyhedron.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Returns the evaluation of shape function i at the supplied point.*/
double PolyhedronMappingFE_PWL::ShapeValue(const int i, const chi_mesh::Vector3& xyz)
{
  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      //Map xyz to xi_eta_zeta
      const auto& p0 = grid->vertices[face_data[f].sides[s].v_index[0]];
      chi_mesh::Vector3 xyz_ref = xyz - p0;

      chi_mesh::Vector3 xi_eta_zeta   = face_data[f].sides[s].Jinv * xyz_ref;

      double xi  = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta= xi_eta_zeta.z;


      //Determine if inside tet
      if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and (zeta>=-1.0e-12) and
          ((xi + eta + zeta)<=(1.0+1.0e-12)))
      {
        double Ni = 0.0;
        double Nf = 0.0;
        double Nc = alphac*zeta;

        if (node_side_maps[i].face_map[f].side_map[s].part_of_face)
        {
          if (node_side_maps[i].face_map[f].side_map[s].index == 0)
          {
            Ni = 1-xi-eta-zeta;
          }
          if (node_side_maps[i].face_map[f].side_map[s].index == 2)
          {
            Ni = eta;
          }

          Nf = face_betaf[f]*xi;
        }

        return Ni + Nf + Nc;
      }
    }
  }
  return 0.0;
}

//###################################################################
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void PolyhedronMappingFE_PWL::ShapeValues(const chi_mesh::Vector3& xyz,
                                          std::vector<double>& shape_values)
{
  shape_values.resize(num_nodes, 0.0);
  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      auto& side_fe_info = face_data[f].sides[s];
      //Map xyz to xi_eta_zeta
      const auto& p0 = grid->vertices[side_fe_info.v_index[0]];
      chi_mesh::Vector3 xi_eta_zeta   = side_fe_info.Jinv * (xyz - p0);

      double xi  = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta= xi_eta_zeta.z;


      //Determine if inside tet
      if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and (zeta>=-1.0e-12) and
          ((xi + eta + zeta)<=(1.0+1.0e-12)))
      {
        for (int i=0; i < num_nodes; i++)
        {
          auto side_map = node_side_maps[i].face_map[f].side_map[s];

          double Ni = 0.0;
          double Nf = 0.0;
          double Nc = alphac*zeta;

          if (side_map.part_of_face)
          {
            if      (side_map.index == 0) Ni = 1-xi-eta-zeta;
            else if (side_map.index == 2) Ni = eta;

            Nf = face_betaf[f]*xi;
          }

          shape_values[i] = Ni + Nf + Nc;
        }//for dof
        return;
      }//if in tet
    }//for side
  }//for face
}

//###################################################################
/**Returns the evaluation of grad-shape function i at the supplied point.*/
chi_mesh::Vector3 PolyhedronMappingFE_PWL::GradShapeValue(const int i,
                                                          const chi_mesh::Vector3& xyz)
{
  chi_mesh::Vector3 grad,gradr;
  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      //Map xyz to xi_eta_zeta
      const auto& p0 = grid->vertices[face_data[f].sides[s].v_index[0]];
      chi_mesh::Vector3 xyz_ref = xyz - p0;

      chi_mesh::Vector3 xi_eta_zeta = face_data[f].sides[s].Jinv * xyz_ref;

      double xi  = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta= xi_eta_zeta.z;

      //Determine if inside tet
      if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and (zeta>=-1.0e-12) and
          ((xi + eta + zeta)<=(1.0+1.0e-12)))
      {
        chi_mesh::Vector3 grad_i;
        chi_mesh::Vector3 grad_f;
        chi_mesh::Vector3 grad_c;

        if (node_side_maps[i].face_map[f].side_map[s].part_of_face)
        {
          if (node_side_maps[i].face_map[f].side_map[s].index == 0)
          {
            grad_i.x =-1.0;
            grad_i.y =-1.0;
            grad_i.z =-1.0;
          }
          if (node_side_maps[i].face_map[f].side_map[s].index == 2)
          {
            grad_i.x = 0.0;
            grad_i.y = 1.0;
            grad_i.z = 0.0;
          }

          grad_f.x = face_betaf[f]*1.0;
          grad_f.y = face_betaf[f]*0.0;
          grad_f.z = face_betaf[f]*0.0;
        }

        grad_c.x = alphac*0.0;
        grad_c.y = alphac*0.0;
        grad_c.z = alphac*1.0;

        grad = (grad_i+grad_f+grad_c);
        grad = face_data[f].sides[s].JTinv * grad;


        return grad;
      }
    }
  }
  return gradr;
}

//###################################################################
/**Populates gradshape_values with the value of each shape function's
 * gradient evaluated at the supplied point.*/
void PolyhedronMappingFE_PWL::GradShapeValues(
  const chi_mesh::Vector3 &xyz,
  std::vector<chi_mesh::Vector3> &gradshape_values)
{
  gradshape_values.clear();
  for (int i=0; i < num_nodes; ++i)
    gradshape_values.emplace_back(GradShapeValue(i,xyz));
}