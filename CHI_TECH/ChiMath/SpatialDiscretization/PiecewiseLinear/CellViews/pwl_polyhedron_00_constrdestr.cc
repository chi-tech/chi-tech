#include "pwl_polyhedron.h"

//###################################################################
/**Constructor for the Piecewise Linear Polyhedron cell finite elment
 * view.
 *
 * */
PolyhedronFEView::PolyhedronFEView(chi_mesh::CellPolyhedron *polyh_cell,
                                   chi_mesh::MeshContinuum *vol_continuum,
                                   SpatialDiscretization_PWL *discretization):
                                   CellFEView(polyh_cell->v_indices.size())
{
  precomputed = false;
  grid = vol_continuum;
  //=========================================== Create quadrature points
  if (discretization != nullptr)
  {
    quadratures.push_back(discretization->tet_quad_deg1);
    quadratures.push_back(discretization->tet_quad_deg3);
    quadratures.push_back(discretization->tet_quad_deg3_surface);
  } else
  {
    fprintf(stderr,"ERROR!: Quadrature points not set for "
                   "pwl_polyhedron.");
    exit(EXIT_FAILURE);
  }


  //=========================================== Assign cell centre
  vcc = polyh_cell->centroid;
  alphac = 1.0/polyh_cell->v_indices.size();

  //=========================================== For each face
  faces.reserve(polyh_cell->faces.size());
  for (int f=0; f<polyh_cell->faces.size(); f++)
  {
    chi_mesh::PolyFace* face = polyh_cell->faces[f];
    FEface_data* face_data = new FEface_data;
    face_betaf.push_back(1.0/face->v_indices.size());

    chi_mesh::Vertex vfc = face->face_centroid;

    face_data->vfc = vfc;

    //==================================== For each edge
    face_data->sides.reserve(face->edges.size());
    for (int e=0; e<face->edges.size(); e++)
    {
      FEside_data3d* side_data = new FEside_data3d;

      //============================= Assign vertices of tetrahedron
      int v0index = face->edges[e][0];
      int v1index = face->edges[e][1];
      side_data->v_index = new int[2];
      side_data->v_index[0] = v0index;
      side_data->v_index[1] = v1index;

      chi_mesh::Vertex v0 = *vol_continuum->nodes[v0index];
      chi_mesh::Vertex v1 = vfc;
      chi_mesh::Vertex v2 = *vol_continuum->nodes[v1index];
      chi_mesh::Vertex v3 = vcc;

      side_data->sc = (v0+v1+v2+v3)/4.0;

      //============================= Compute vectors
      chi_mesh::Vector v01 = v1 - v0;
      chi_mesh::Vector v02 = v2 - v0;
      chi_mesh::Vector v03 = v3 - v0;

      //============================= Compute determinant of surface jacobian
      // First we compute the rotation matrix which will rotate
      // any vector in natural coordinates to the same reference
      // frame as the current face.
      chi_mesh::Vector normal  = face->geometric_normal*-1;
      chi_mesh::Vector tangent = v02.Cross(normal);
      tangent = tangent/tangent.Norm();
      chi_mesh::Vector binorm  = v02/v02.Norm();

      chi_mesh::Matrix3x3 R;
      R.SetColJVec(0,tangent);
      R.SetColJVec(1,binorm);
      R.SetColJVec(2,normal);

      // Now we compute the inverse of this matrix which
      // will allow us to rotate any vector in the same reference
      // frame as the face, to natural coordinates
      chi_mesh::Matrix3x3 Rinv = R.Inverse();

      // Compute v01 and v02 rotated to natural coordinates
      // A test to see if this is done correctly would be to
      // check if fabs(v01N.z) < epsilon and fabs(v02N.z) < epsilon
      chi_mesh::Vector v01N = Rinv*v01;
      chi_mesh::Vector v02N = Rinv*v02;
//        if ( (std::fabs(v01N.z) > 1.0e-5) || (std::fabs(v01N.z) > 1.0e-5) )
//        {
//          printf("ERROR! v01\n");
//        }
      side_data->detJ_surf = v01N.x*v02N.y - v01N.y*v02N.x;

      //============================= Compute Jacobian
      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0,v01);
      J.SetColJVec(1,v02);
      J.SetColJVec(2,v03);

      side_data->J = J;


      //============================= Compute determinant of jacobian
      side_data->detJ = J.Det();

      //============================= Compute inverse Jacobian elements
      chi_mesh::Matrix3x3 JT   = J.Transpose();
      chi_mesh::Matrix3x3 Jinv = J.Inverse();
      chi_mesh::Matrix3x3 JTinv= JT.Inverse();

      side_data->Jinv  = Jinv;
      side_data->JTinv = JTinv;

      side_data->qp_data.reserve(dofs);
      face_data->sides.push_back(side_data);
    }//for each edge

    faces.push_back(face_data);
  }//for each face


  //================================================ Compute Node-Face-Side
  //                                                 mapping
  // This section determines the scope of dof_i on
  // each side (tet) of the cell. If dof_i is on
  // either of the primary tet nodes, it is given
  // index 0 or 1 (-1 otherwise). If the index is
  // -1 the corresponding shapefunction will be
  // ignored when using Ni. The flag "part_of_face"
  // is set when dof_i is part of the same face to
  // which the side belongs and consequently allows
  // the determination of Nf. Nc is always evaluated
  // so no mapping is needed.
  for (int i=0; i<dofs; i++)
  {
    FEnodeMap* newNodeMap = new FEnodeMap;
    for (int f=0; f<faces.size(); f++)
    {
      FEnodeFaceMap* newFaceMap = new FEnodeFaceMap;
      for (int s=0; s<faces[f]->sides.size(); s++)
      {
        FEnodeSideMap* newSideMap = new FEnodeSideMap;
        newSideMap->part_of_face = false;
        int s0 = faces[f]->sides[s]->v_index[0];
        int s1 = faces[f]->sides[s]->v_index[1];
        if      (polyh_cell->v_indices[i] == s0)
        {
          newSideMap->index = 0;
          newSideMap->part_of_face = true;
        }
        else if (polyh_cell->v_indices[i] == s1)
        {
          newSideMap->index = 2;
          newSideMap->part_of_face = true;
        }
        else
        {
          newSideMap->index = -1;
          for (int v=0; v<polyh_cell->faces[f]->v_indices.size(); v++)
          {
            if (polyh_cell->v_indices[i] ==
                polyh_cell->faces[f]->v_indices[v])
            {
              newSideMap->part_of_face = true;
              break;
            }
          }

        }
        newFaceMap->side_map.push_back(newSideMap);
      }//for s
      newNodeMap->face_map.push_back(newFaceMap);
    }//for f
    node_maps.push_back(newNodeMap);
  }//for i

  //================================================ Compute Face DOF mapping
  // This section just determines a mapping of face dofs
  // to cell dofs. This is pretty simple since we can
  // just loop over each face dof then subsequently
  // loop over cell dofs, if the face dof node index equals
  // the cell dof node index then the mapping is assigned.
  //
  // This mapping is not used by any of the methods in
  // this class but is used by methods requiring the
  // surface integrals of the shape functions.
  for (int f=0; f<polyh_cell->faces.size(); f++)
  {
    FEFaceDOFMapping* face_dof_mapping = new FEFaceDOFMapping;

    for (int fi=0; fi<polyh_cell->faces[f]->v_indices.size(); fi++)
    {
      int mapping = -1;
      for (int ci=0; ci<polyh_cell->v_indices.size(); ci++)
      {
        if (polyh_cell->faces[f]->v_indices[fi] ==
            polyh_cell->v_indices[ci])
        {
          mapping = ci;
          break;
        }
      }//for cell i
      if (mapping<0)
      {
        chi_log.Log(LOG_ALLERROR) << "Unknown face mapping encountered. "
                                     "pwl_polyhedron.h";
        exit(EXIT_FAILURE);
      }
      face_dof_mapping->cell_dof.push_back(mapping);
    }//for face i

    face_dof_mappings.push_back(face_dof_mapping);
  }
}