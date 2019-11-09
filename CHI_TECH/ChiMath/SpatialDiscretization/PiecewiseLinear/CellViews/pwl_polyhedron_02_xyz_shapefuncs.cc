#include "pwl_polyhedron.h"

double PolyhedronFEView::Shape_xyz(int i, chi_mesh::Vector& xyz)
{
  for (int f=0; f<faces.size(); f++)
  {
    for (int s=0; s<faces[f]->sides.size(); s++)
    {
      //Map xyz to xi_eta_zeta
      chi_mesh::Vector p0 = *grid->nodes[faces[f]->sides[s]->v_index[0]];
      chi_mesh::Vector xyz_ref = xyz - p0;

      chi_mesh::Vector xi_eta_zeta   = faces[f]->sides[s]->Jinv*xyz_ref;

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

        if (node_maps[i]->face_map[f]->side_map[s]->part_of_face)
        {
          if (node_maps[i]->face_map[f]->side_map[s]->index == 0)
          {
            Ni = 1-xi-eta-zeta;
          }
          if (node_maps[i]->face_map[f]->side_map[s]->index == 2)
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

chi_mesh::Vector PolyhedronFEView::GradShape_xyz(int i, chi_mesh::Vector xyz)
{
  chi_mesh::Vector grad,gradr;
  for (int f=0; f<faces.size(); f++)
  {
    for (int s=0; s<faces[f]->sides.size(); s++)
    {
      //Map xyz to xi_eta_zeta
      chi_mesh::Vector p0 = *grid->nodes[faces[f]->sides[s]->v_index[0]];
      chi_mesh::Vector xyz_ref = xyz - p0;

      chi_mesh::Vector xi_eta_zeta = faces[f]->sides[s]->Jinv*xyz_ref;

      double xi  = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta= xi_eta_zeta.z;

      //Determine if inside tet
      if ((xi>=-1.0e-12) and (eta>=-1.0e-12) and (zeta>=-1.0e-12) and
          ((xi + eta + zeta)<=(1.0+1.0e-12)))
      {
        chi_mesh::Vector grad_i;
        chi_mesh::Vector grad_f;
        chi_mesh::Vector grad_c;

        if (node_maps[i]->face_map[f]->side_map[s]->part_of_face)
        {
          if (node_maps[i]->face_map[f]->side_map[s]->index == 0)
          {
            grad_i.x =-1.0;
            grad_i.y =-1.0;
            grad_i.z =-1.0;
          }
          if (node_maps[i]->face_map[f]->side_map[s]->index == 2)
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
        grad = faces[f]->sides[s]->JTinv*grad;


        return grad;
      }
    }
  }
  return gradr;
}