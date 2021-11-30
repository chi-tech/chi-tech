#include "pwl_polygon.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/** Constructor.*/
PolygonMappingFE_PWL::
  PolygonMappingFE_PWL(
    const chi_mesh::Cell& poly_cell,
    const std::shared_ptr<chi_mesh::MeshContinuum>& ref_grid,
    const chi_math::QuadratureTriangle& volume_quadrature,
    const chi_math::QuadratureLine&     surface_quadrature) :
  CellMappingFE_PWL(poly_cell.vertex_ids.size(), ref_grid),
  volume_quadrature(volume_quadrature),
  surface_quadrature(surface_quadrature)
{
  num_of_subtris = poly_cell.faces.size();
  beta = 1.0/num_of_subtris;

  //=========================================== Get raw vertices
  vc = poly_cell.centroid;

  //=========================================== Calculate legs and determinants
  for (int side=0;side<num_of_subtris;side++)
  {
    const chi_mesh::CellFace& face = poly_cell.faces[side];

    const auto& v0 = ref_grid->vertices[face.vertex_ids[0]];
    const auto& v1 = ref_grid->vertices[face.vertex_ids[1]];
    chi_mesh::Vertex v2 = vc;

    chi_mesh::Vector3 sidev01 = v1 - v0;
    chi_mesh::Vector3 sidev02 = v2 - v0;

    double sidedetJ = ((sidev01.x)*(sidev02.y) - (sidev02.x)*(sidev01.y));

    FEside_data2d triangle_data;
    triangle_data.detJ = sidedetJ;
    triangle_data.detJ_surf = sidev01.Norm();

    triangle_data.v_index[0] = face.vertex_ids[0];
    triangle_data.v_index[1] = face.vertex_ids[1];

    triangle_data.v0 = v0;

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

    //Set face normal
    triangle_data.normal = face.normal;

    sides.push_back(triangle_data);
  }

  //=========================================== Compute node to side mapping
  for (int v=0; v<poly_cell.vertex_ids.size(); v++)
  {
    int vindex = poly_cell.vertex_ids[v];
    std::vector<int> side_mapping(num_of_subtris);
    for (int side=0;side<num_of_subtris;side++)
    {
      side_mapping[side] = -1;

      const chi_mesh::CellFace& face = poly_cell.faces[side];
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
  face_dof_mappings.resize(poly_cell.faces.size());
  for (int e=0; e<poly_cell.faces.size(); e++)
  {
    face_dof_mappings[e].resize(2);
    for (int fv=0; fv<2; fv++)
    {
      for (int v=0; v<poly_cell.vertex_ids.size(); v++)
      {
        if (poly_cell.faces[e].vertex_ids[fv] == poly_cell.vertex_ids[v])
        {
          face_dof_mappings[e][fv] = v;
          break;
        }
      }//for v
    }//for fv
  }//for e
}