#include "chi_ffinter_line.h"

#include "../../Cell/cell_slab.h"
#include "../../Cell/cell_polygon.h"
#include "../../Cell/cell_polyhedron.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Initializes the data structures necessary for interpolation. This is
 * independent of the physics and hence is a routine on its own.
 *
 * The first step of this initialization is to determine which cells
 * are intersected by this line. For polyhedrons this is evaluated
 * face side-by-side.
 *
 * The second step is to find the cell associated with each point in
 * in the line.
 *
 * Third step is to upload node indices for each-cell-point pair so that
 * the value can be interpolated.*/
void chi_mesh::FieldFunctionInterpolationLine::
Initialize()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing line interpolator2.";
  //================================================== Check grid available
  if (field_functions.size() == 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unassigned field function in line field function interpolator.";
    exit(EXIT_FAILURE);
  }
  grid_view = nullptr;

  //================================================== Create points;
  chi_mesh::Vector vif = pf-pi;
  delta_d = vif.Norm()/(number_of_points-1);

  vif = vif/vif.Norm();

  interpolation_points.push_back(pi);
  for (int k=1; k<(number_of_points); k++)
  {
    interpolation_points.push_back(pi + vif*delta_d*k);
  }

  //====================================================== Loop over contexts
  for (int ff=0; ff<field_functions.size(); ff++)
  {
    auto ff_context = new FieldFunctionContext;
    ff_contexts.push_back(ff_context);

    ff_context->ref_ff = field_functions[ff];

    if (grid_view != field_functions[ff]->grid)
    {
      grid_view = field_functions[ff]->grid;
      interpolation_points_ass_cell.resize(number_of_points,-1);
      interpolation_points_ass_cell.assign(number_of_points,-1);
      cfem_local_nodes_needed_unmapped.clear();
      pwld_local_nodes_needed_unmapped.clear();
      pwld_local_cells_needed_unmapped.clear();

      //================================================== Find a home for each
      //                                                   point
      size_t num_local_cells = grid_view->local_cell_glob_indices.size();
      for (int ic=0; ic<num_local_cells; ic++)
      {
        int cell_glob_index = grid_view->local_cell_glob_indices[ic];
        auto cell = grid_view->cells[cell_glob_index];

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
        if (cell->Type() == chi_mesh::CellType::SLAB)
        {
          chi_mesh::CellSlab* slab_cell = (chi_mesh::CellSlab*)cell;



          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            int v0_i = slab_cell->v_indices[0];
            int v1_i = slab_cell->v_indices[1];

            chi_mesh::Vector v0 = *grid_view->nodes[v0_i];
            chi_mesh::Vector v1 = *grid_view->nodes[v1_i];

            chi_mesh::Vector v01 = v1 - v0;
            chi_mesh::Vector v0p = interpolation_points[p]-v0;

            double norm = v0p.Norm();
            double projection = v01.Dot(v0p)/v01.Norm();

            if ((v0p.Dot(v01)<0.0) or (projection>v01.Norm()))
            {
              is_inside = false;
            }

            if (is_inside and interpolation_points_ass_cell[p]<0)
            {
              interpolation_points_ass_cell[p] = cell_glob_index;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p;
            }

          }//for each point
        }//if slab

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
        else if (cell->Type() == chi_mesh::CellType::POLYGON)
        {
          chi_mesh::CellPolygon* poly_cell = (chi_mesh::CellPolygon*)cell;

          size_t num_edges = poly_cell->edges.size();

          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            //Form a plane for each edge
            for (int e=0; e<num_edges; e++)
            {
              chi_mesh::Vector nref(0.0,0.0,1.0);
              int v0_i = poly_cell->edges[e][0];
              int v1_i = poly_cell->edges[e][1];

              chi_mesh::Vector v0 = *grid_view->nodes[v0_i];
              chi_mesh::Vector v1 = *grid_view->nodes[v1_i];

              chi_mesh::Vector v01 = v1 - v0;
              chi_mesh::Vector   n = v01.Cross(nref);
              n = n/n.Norm();

              chi_mesh::Vector v0p = interpolation_points[p]-v0;
              v0p=v0p/v0p.Norm();

              if (n.Dot(v0p)>0.0)
              {
                is_inside = false;
                break;
              }
            }//for edge

            if (is_inside)
            {
              interpolation_points_ass_cell[p] = cell_glob_index;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p;
            }

          }//for each point
        }//if polygon cell

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
        else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
        {
          chi_mesh::CellPolyhedron* polyh_cell = (chi_mesh::CellPolyhedron*)cell;

          size_t num_faces = polyh_cell->faces.size();

          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            for (int f=0; f<num_faces; f++)
            {
              size_t num_edges = polyh_cell->faces[f]->edges.size();
              //Form a plane for each side
              for (int e=0; e<num_edges; e++)
              {
                int v0_i = polyh_cell->faces[f]->edges[e][0];

                chi_mesh::Vector v0 = *grid_view->nodes[v0_i];
                chi_mesh::Vector n  = polyh_cell->faces[f]->geometric_normal;

                chi_mesh::Vector v0p = interpolation_points[p]-v0;
                v0p=v0p/v0p.Norm();

                if (n.Dot(v0p)>0.0)
                {
                  is_inside = false;
                  break;
                }
              }//for e
              if (!is_inside) break;
            }//for f

            if (is_inside)
            {
              interpolation_points_ass_cell[p] = cell_glob_index;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p;
            }

          }//for each point

        }//if polyh
      }//for local cell

      //================================================== Upload node indices that
      //                                                   need mapping
      for (int c=0; c<interpolation_points_ass_cell.size(); c++)
      {
        if (interpolation_points_ass_cell[c] < 0) continue;

        int cell_glob_index = interpolation_points_ass_cell[c];
        auto cell = grid_view->cells[cell_glob_index];

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
        if (cell->Type() == chi_mesh::CellType::SLAB)
        {
          chi_mesh::CellSlab* slab_cell = (chi_mesh::CellSlab*)cell;

          for (int i=0; i<2; i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(slab_cell->v_indices[i]);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
          }

        }//if poly

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
        if (cell->Type() == chi_mesh::CellType::POLYGON)
        {
          chi_mesh::CellPolygon* poly_cell = (chi_mesh::CellPolygon*)cell;

          for (int i=0; i<poly_cell->v_indices.size(); i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(poly_cell->v_indices[i]);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
          }

        }//if poly

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
        else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
        {
          chi_mesh::CellPolyhedron* polyh_cell = (chi_mesh::CellPolyhedron*)cell;

          for (int i=0; i<polyh_cell->v_indices.size(); i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(polyh_cell->v_indices[i]);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_glob_index);
          }

        }//if polyh
      }//for associated cells
    }//if unique grid


    //Copies the latest developed references to the specific
    //field function context
    std::copy(cfem_local_nodes_needed_unmapped.begin(),
              cfem_local_nodes_needed_unmapped.end(),
              std::back_inserter(
              ff_context->cfem_local_nodes_needed_unmapped));
    std::copy(pwld_local_nodes_needed_unmapped.begin(),
              pwld_local_nodes_needed_unmapped.end(),
              std::back_inserter(
                ff_context->pwld_local_nodes_needed_unmapped));
    std::copy(pwld_local_cells_needed_unmapped.begin(),
              pwld_local_cells_needed_unmapped.end(),
              std::back_inserter(
                ff_context->pwld_local_cells_needed_unmapped));
    std::copy(interpolation_points_ass_cell.begin(),
              interpolation_points_ass_cell.end(),
              std::back_inserter(
                ff_context->interpolation_points_ass_cell));
  }//for ff



  chi_log.Log(LOG_0VERBOSE_1) << "Finished initializing interpolator.";
}