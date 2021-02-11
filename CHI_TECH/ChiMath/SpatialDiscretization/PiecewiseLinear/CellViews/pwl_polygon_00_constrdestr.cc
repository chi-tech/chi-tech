#include "pwl_polygon.h"

//###################################################################
/** Constructor.*/
PolygonPWLFEValues::PolygonPWLFEValues(chi_mesh::CellPolygon* poly_cell,
                                       chi_mesh::MeshContinuumPtr vol_continuum,
                                       SpatialDiscretization_PWL *discretization) :
  CellPWLFEValues(poly_cell->vertex_ids.size()),
  default_volume_quadrature(discretization->tri_quad_order_second),
  default_surface_quadrature(discretization->line_quad_order_second)
{
  precomputed = false;
  grid = vol_continuum;
  num_of_subtris = poly_cell->faces.size();
  beta = 1.0/num_of_subtris;

  //=========================================== Get raw vertices
  vc = poly_cell->centroid;

  //=========================================== Calculate legs and determinants
  for (int side=0;side<num_of_subtris;side++)
  {
    chi_mesh::CellFace& face = poly_cell->faces[side];

    chi_mesh::Vertex v0 = *vol_continuum->vertices[face.vertex_ids[0]];
    chi_mesh::Vertex v1 = *vol_continuum->vertices[face.vertex_ids[1]];
    chi_mesh::Vertex v2 = vc;

    chi_mesh::Vector3 sidev01 = v1 - v0;
    chi_mesh::Vector3 sidev02 = v2 - v0;

    double sidedetJ = ((sidev01.x)*(sidev02.y) - (sidev02.x)*(sidev01.y));
    detJ.push_back(sidedetJ);

    FEside_data2d triangle_data;
    triangle_data.detJ = sidedetJ;
    triangle_data.detJ_surf = sidev01.Norm();

    triangle_data.v_index[0] = face.vertex_ids[0];
    triangle_data.v_index[1] = face.vertex_ids[1];

    //Set Jacobian
    triangle_data.J.SetIJ(0,0,sidev01.x);
    triangle_data.J.SetIJ(1,0,sidev01.y);
    triangle_data.J.SetIJ(0,1,sidev02.x);
    triangle_data.J.SetIJ(1,1,sidev02.y);
    triangle_data.J.SetIJ(2,2,0.0);

    //Set Jacobian inverse
    triangle_data.Jinv.SetIJ(0,0,sidev02.y/sidedetJ);
    triangle_data.Jinv.SetIJ(1,0,-sidev01.y/sidedetJ);
    triangle_data.Jinv.SetIJ(0,1,-sidev02.x/sidedetJ);
    triangle_data.Jinv.SetIJ(1,1,sidev01.x/sidedetJ);
    triangle_data.Jinv.SetIJ(2,2,0.0);

    //Set Jacobian-Transpose inverse
    triangle_data.JTinv.SetIJ(0,0, sidev02.y/sidedetJ);
    triangle_data.JTinv.SetIJ(1,0,-sidev02.x/sidedetJ);
    triangle_data.JTinv.SetIJ(0,1,-sidev01.y/sidedetJ);
    triangle_data.JTinv.SetIJ(1,1, sidev01.x/sidedetJ);
    triangle_data.JTinv.SetIJ(2,2,0.0);

    sides.push_back(triangle_data);
  }

  //=========================================== Compute node to side mapping
  for (int v=0; v<poly_cell->vertex_ids.size(); v++)
  {
    int vindex = poly_cell->vertex_ids[v];
    int* side_mapping = new int[num_of_subtris];
    for (int side=0;side<num_of_subtris;side++)
    {
      side_mapping[side] = -1;

      chi_mesh::CellFace& face = poly_cell->faces[side];
      if (face.vertex_ids[0] == vindex)
      {
        side_mapping[side] = 0;
      }
      if (face.vertex_ids[1] == vindex)
      {
        side_mapping[side] = 1;
      }
    }
    node_to_side_map.push_back(side_mapping);
  }

  //============================================= Compute edge dof mappings
  face_dof_mappings.resize(poly_cell->faces.size());
  for (int e=0; e<poly_cell->faces.size(); e++)
  {
    face_dof_mappings[e].resize(2);
    for (int fv=0; fv<2; fv++)
    {
      for (int v=0; v<poly_cell->vertex_ids.size(); v++)
      {
        if (poly_cell->faces[e].vertex_ids[fv] == poly_cell->vertex_ids[v])
        {
          face_dof_mappings[e][fv] = v;
          break;
        }
      }
    }
  }
}