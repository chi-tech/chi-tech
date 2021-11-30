#include "chi_ffinter_slice.h"
#include "ChiMesh/Cell/cell.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include <chi_log.h>

extern ChiLog& chi_log;

/**Initializes the data structures necessary for interpolation. This is
 * independent of the physics and hence is a routine on its own.
 *
 * The first step of this initialization is to determine which cells
 * are intersected by this plane. For polyhedrons this is evaluated
 * tet-by-tet.
 *
 * The second step is find where face-edges are intersected. This will
 * effectively create intersection polygons.*/
void chi_mesh::FieldFunctionInterpolationSlice::
  Initialize()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing slice interpolator.";
  //================================================== Check grid available
  if (field_functions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
    << "Unassigned field function in slice field function interpolator.";
    exit(EXIT_FAILURE);
  } else
  {
    this->grid_view = field_functions[0]->grid;
  }

  //================================================== Find cells intersecting plane
  intersecting_cell_indices.clear();

  for (const auto& cell : grid_view->local_cells)
  {
    auto cell_local_index = cell.local_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      chi_log.Log(LOG_0)
        << "FieldFunctionInterpolationSlice does not support 1D cells.";
      exit(EXIT_FAILURE);
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      intersecting_cell_indices.push_back(cell_local_index);
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      bool intersects = false;

      size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; f++)
      {
        const auto& face = cell.faces[f];

        size_t num_edges = face.vertex_ids.size();
        for (size_t e=0; e<num_edges; e++)
        {
          size_t ep1 = (e < (num_edges-1))? e+1 : 0;
          uint64_t v0_i = face.vertex_ids[e];
          uint64_t v1_i = face.vertex_ids[ep1];

          std::vector<chi_mesh::Vector3> tet_points;

          tet_points.push_back(grid_view->vertices[v0_i]);
          tet_points.push_back(grid_view->vertices[v1_i]);
          tet_points.push_back(cell.faces[f].centroid);
          tet_points.push_back(cell.centroid);

          if (CheckPlaneTetIntersect(this->normal,this->point,tet_points))
          {
            intersecting_cell_indices.push_back(cell_local_index);
            intersects = true;
            break;
          }
        }//for e
        if (intersects) break;
      }//for f
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type in call to Slice Initialize.";
      exit(EXIT_FAILURE);
    }
  }//for local cell

  //================================================== Computing cell intersections
  size_t num_cut_cells = intersecting_cell_indices.size();
  for (int cc=0; cc<num_cut_cells; cc++)
  {
    uint64_t cell_local_index = intersecting_cell_indices[cc];

    const auto& cell = grid_view->local_cells[cell_local_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      chi_log.Log(LOG_0)
        << "FieldFunctionInterpolationSlice does not support 1D cells.";
      exit(EXIT_FAILURE);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      //========================================= Initialize cell intersection
      //                                          data structure
      FFICellIntersection cell_isds;
      cell_isds.cell_local_index = cell_local_index;

      //========================================= Loop over vertices
      for (int v=0; v<cell.vertex_ids.size(); v++)
      {
        FFIFaceEdgeIntersection face_isds;

        int v0gi = cell.vertex_ids[v];

        face_isds.v0_g_index = v0gi;
        face_isds.v1_g_index = v0gi;

        face_isds.v0_dofindex_cell = v;
        face_isds.v1_dofindex_cell = v;

        face_isds.weights = std::pair<double,double>(0.5,0.5);
        face_isds.point   = grid_view->vertices[v0gi];

        cell_isds.intersections.push_back(face_isds);
      }

      //========================================= Set intersection center
      cell_isds.intersection_centre = cell.centroid;

      //========================================= Set straight 2D center
      // This is normally transformed for the 3D case
      cell_isds.intersection_2d_centre = cell_isds.intersection_centre;

      //========================================= Same for 2D points
      size_t num_points = cell_isds.intersections.size();
      for (size_t p=0; p<num_points; p++)
      {
        chi_mesh::Vector3 vref = cell_isds.intersections[p].point - this->point;

        cell_isds.intersections[p].point2d = vref;

        cfem_local_nodes_needed_unmapped.push_back(cell_isds.intersections[p].v0_dofindex_cell);
        cfem_local_nodes_needed_unmapped.push_back(cell_isds.intersections[p].v1_dofindex_cell);
        cfem_local_cells_needed_unmapped.push_back(cell_local_index);
        cfem_local_cells_needed_unmapped.push_back(cell_local_index);

        pwld_local_nodes_needed_unmapped.push_back(cell_isds.intersections[p].v0_dofindex_cell);
        pwld_local_nodes_needed_unmapped.push_back(cell_isds.intersections[p].v1_dofindex_cell);
        pwld_local_cells_needed_unmapped.push_back(cell_local_index);
        pwld_local_cells_needed_unmapped.push_back(cell_local_index);
      }

      cell_intersections.push_back(cell_isds);
    }//polygon
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      //========================================= Initialize cell intersection
      //                                          data structure
      FFICellIntersection cell_isds;
      cell_isds.cell_local_index = cell_local_index;

      //========================================= Loop over faces
      size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; f++)
      {
        const auto& face = cell.faces[f];
        //================================== Loop over edges
        size_t num_edges = face.vertex_ids.size();
        for (size_t e=0; e<num_edges; e++)
        {
          size_t ep1 = (e < (num_edges-1))? e+1 : 0;
          uint64_t v0gi = face.vertex_ids[e  ]; //global index v0
          uint64_t v1gi = face.vertex_ids[ep1]; //global index v1

          const auto& v0 = grid_view->vertices[v0gi];
          const auto& v1 = grid_view->vertices[v1gi];

          chi_mesh::Vertex interstion_point;            //Placeholder
          std::pair<double,double> weights;

          //=========================== Check if intersects plane
          if (CheckPlaneLineIntersect(this->normal,this->point,
                                      v0,v1,interstion_point,
                                      &weights))
          {
            //==================== Check for duplicate
            bool duplicate_found = false;
            for (int is=0; is<cell_isds.intersections.size(); is++)
            {
              FFIFaceEdgeIntersection& existing_face_isds =
                cell_isds.intersections[is];

              double dif = (existing_face_isds.point -
                            interstion_point).NormSquare();
              if (dif < 1.0e-6)
              {
                duplicate_found = true;
                break;
              }
            }

            //==================== No duplicate
            if (!duplicate_found)
            {
              FFIFaceEdgeIntersection face_isds;

              //Find vertex 0 dof index
              face_isds.v0_g_index = v0gi;
              size_t num_dofs = cell.vertex_ids.size();
              for (int dof=0; dof<num_dofs; dof++)
              {
                if (cell.vertex_ids[dof] == v0gi)
                {
                  face_isds.v0_dofindex_cell = dof; break;
                }
              }

              //Find vertex 1 dof index
              face_isds.v1_g_index = v1gi;
              for (int dof=0; dof<num_dofs; dof++)
              {
                if (cell.vertex_ids[dof] == v1gi)
                {
                  face_isds.v1_dofindex_cell = dof; break;
                }
              }

              face_isds.weights = weights;
              face_isds.point = interstion_point;
              cell_isds.intersections.push_back(face_isds);
            }
          }//if intersecting
        }//for edge
      }//for face

      //====================================

      //==================================== Computing intersection centre
      size_t num_points = cell_isds.intersections.size();
      if (num_points>0)
      {
        for (int p=0; p<num_points; p++)
        {
          cell_isds.intersection_centre =
            cell_isds.intersection_centre +
            cell_isds.intersections[p].point;
        }
        cell_isds.intersection_centre =
          cell_isds.intersection_centre/num_points;
      }
      else
      {
        chi_log.Log(LOG_ALLWARNING) << "No face intersections encountered "
                                       "for a cell that is indicated as being "
                                       "intersected. Slice FF interp.";
      }

      //==================================== Computing 2D transforms
      chi_mesh::Vector3 vref = cell_isds.intersection_centre - this->point;

      cell_isds.intersection_2d_centre.x = vref.Dot(tangent);
      cell_isds.intersection_2d_centre.y = vref.Dot(binorm);
      cell_isds.intersection_2d_centre.z = vref.Dot(normal);

      //==================================== Points
      std::vector<FFIFaceEdgeIntersection> unsorted_points;
      for (int p=0; p<num_points; p++)
      {
        chi_mesh::Vector3 vref = cell_isds.intersections[p].point - this->point;

        cell_isds.intersections[p].point2d.x = vref.Dot(tangent);
        cell_isds.intersections[p].point2d.y = vref.Dot(binorm);
        cell_isds.intersections[p].point2d.z = vref.Dot(normal);

        unsorted_points.push_back(cell_isds.intersections[p]);
      }
      cell_isds.intersections.clear();

      //==================================== Sort points clockwise
      //The first point is retrieved from the unused stack.
      //Subsequent points are only added if they form a
      //convex line wrt the right hand rule.
      cell_isds.intersections.push_back(unsorted_points[0]);
      cfem_local_nodes_needed_unmapped.push_back(unsorted_points[0].v0_dofindex_cell);
      cfem_local_nodes_needed_unmapped.push_back(unsorted_points[0].v1_dofindex_cell);
      cfem_local_cells_needed_unmapped.push_back(cell_local_index);
      cfem_local_cells_needed_unmapped.push_back(cell_local_index);

      pwld_local_nodes_needed_unmapped.push_back(unsorted_points[0].v0_dofindex_cell);
      pwld_local_nodes_needed_unmapped.push_back(unsorted_points[0].v1_dofindex_cell);
      pwld_local_cells_needed_unmapped.push_back(cell_local_index);
      pwld_local_cells_needed_unmapped.push_back(cell_local_index);
      unsorted_points.erase(unsorted_points.begin());

      while (unsorted_points.size()>0)
      {
        for (int p=0; p<unsorted_points.size(); p++)
        {
          chi_mesh::Vector3 v1 = unsorted_points[p].point2d -
                                 cell_isds.intersections.back().point2d;

          bool illegal_value = false;
          for (int pr=0; pr<unsorted_points.size(); pr++)
          {
            if (pr!=p)
            {
              chi_mesh::Vector3 vr = unsorted_points[pr].point2d -
                                     unsorted_points[p].point2d;

              if (vr.Cross(v1).z < 0.0)
              {
                illegal_value = true;
                break;
              }
            }//if not p
          }//for pr

          if (!illegal_value)
          {
            cell_isds.intersections.push_back(unsorted_points[p]);
            cfem_local_nodes_needed_unmapped.push_back(unsorted_points[p].v0_dofindex_cell);
            cfem_local_nodes_needed_unmapped.push_back(unsorted_points[p].v1_dofindex_cell);
            cfem_local_cells_needed_unmapped.push_back(cell_local_index);
            cfem_local_cells_needed_unmapped.push_back(cell_local_index);

            pwld_local_nodes_needed_unmapped.push_back(unsorted_points[p].v0_dofindex_cell);
            pwld_local_nodes_needed_unmapped.push_back(unsorted_points[p].v1_dofindex_cell);
            pwld_local_cells_needed_unmapped.push_back(cell_local_index);
            pwld_local_cells_needed_unmapped.push_back(cell_local_index);
            unsorted_points.erase(unsorted_points.begin()+p);
            break;
          }

        }//for p
      }
      cell_intersections.push_back(cell_isds);
    }//polyhedron
  }//for intersected cell

  //chi_log.Log(LOG_0) << "Finished initializing interpolator.";
}