#include "chi_ffinter_line.h"

#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
extern ChiLog& chi_log;

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
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing line interpolator.";
  //================================================== Check grid available
  if (field_functions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unassigned field function in line field function interpolator.";
    exit(EXIT_FAILURE);
  }
  grid_view = nullptr;

  //================================================== Create points;
  chi_mesh::Vector3 vif = pf - pi;
  delta_d = vif.Norm()/(number_of_points-1);

  vif = vif/vif.Norm();

  interpolation_points.push_back(pi);
  for (int k=1; k<(number_of_points); k++)
    interpolation_points.push_back(pi + vif*delta_d*k);

  //====================================================== Initialize scratch data
  std::vector<int>               cfem_local_nodes_needed_unmapped;
  std::vector<int>               cfem_local_cells_needed_unmapped;
  std::vector<int>               pwld_local_nodes_needed_unmapped;
  std::vector<int>               pwld_local_cells_needed_unmapped;
  std::vector<uint64_t>          interpolation_points_ass_cell;
  std::vector<bool>              interpolation_points_has_ass_cell;

  //====================================================== Loop over contexts
  size_t num_ff = field_functions.size();
  for (size_t ff=0; ff<num_ff; ff++)
  {
    auto ff_context = new FieldFunctionContext;
    ff_contexts.push_back(ff_context);

    ff_context->ref_ff = field_functions[ff];

    if (grid_view != field_functions[ff]->grid)
    {
      grid_view = field_functions[ff]->grid;
      interpolation_points_ass_cell.resize(number_of_points,0);
      interpolation_points_ass_cell.assign(number_of_points,0);

      interpolation_points_has_ass_cell.resize(number_of_points,false);
      interpolation_points_has_ass_cell.assign(number_of_points,false);

      cfem_local_nodes_needed_unmapped.clear();
      cfem_local_cells_needed_unmapped.clear();
      pwld_local_nodes_needed_unmapped.clear();
      pwld_local_cells_needed_unmapped.clear();

      //================================================== Find a home for each
      //                                                   point
      for (const auto& cell : grid_view->local_cells)
      {
        int cell_local_index = cell.local_id;

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
        if (cell.Type() == chi_mesh::CellType::SLAB)
        {
          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            uint64_t v0_i = cell.vertex_ids[0];
            uint64_t v1_i = cell.vertex_ids[1];

            const auto& v0 = grid_view->vertices[v0_i];
            const auto& v1 = grid_view->vertices[v1_i];

            chi_mesh::Vector3 v01 = v1 - v0;
            chi_mesh::Vector3 v0p = interpolation_points[p] - v0;

            double v01_norm = v01.Norm();
            double projection = v01.Dot(v0p)/v01_norm;

            if ((v0p.Dot(v01)<0.0) or (projection>v01_norm))
            {
              is_inside = false;
            }

            if (is_inside)
            {
              interpolation_points_ass_cell[p] = cell_local_index;
              interpolation_points_has_ass_cell[p] = true;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p
                << " v0[" << v0_i << "]=" << v0.PrintS()
                << " v1[" << v1_i << "]=" << v1.PrintS()
                << " point=" << interpolation_points[p].PrintS();

            }

          }//for each point
        }//if slab

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
        else if (cell.Type() == chi_mesh::CellType::POLYGON)
        {
          size_t num_edges = cell.faces.size();

          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            //Form a plane for each edge
            for (size_t e=0; e<num_edges; e++)
            {
              chi_mesh::Vector3 nref(0.0, 0.0, 1.0);
              uint64_t v0_i = cell.faces[e].vertex_ids[0];
              uint64_t v1_i = cell.faces[e].vertex_ids[1];

              const auto& v0 = grid_view->vertices[v0_i];
              const auto& v1 = grid_view->vertices[v1_i];

              chi_mesh::Vector3 v01 = v1 - v0;
              chi_mesh::Vector3   n = v01.Cross(nref);
              n = n/n.Norm();

              chi_mesh::Vector3 v0p = interpolation_points[p] - v0;
              v0p=v0p/v0p.Norm();

              if (n.Dot(v0p)>0.0)
              {
                is_inside = false;
                break;
              }
            }//for edge

            if (is_inside)
            {
              interpolation_points_ass_cell[p] = cell_local_index;
              interpolation_points_has_ass_cell[p] = true;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p;
            }

          }//for each point
        }//if polygon cell

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
        else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
        {
          size_t num_faces = cell.faces.size();

          for (int p=0; p<number_of_points; p++)
          {
            //Assume each point is inside the cell, now try to disprove it
            bool is_inside = true;
            for (size_t f=0; f<num_faces; f++)
            {
              const auto& face = cell.faces[f];
              size_t num_edges = face.vertex_ids.size();
              //Form a plane for each side
              for (size_t e=0; e<num_edges; e++)
              {
                uint64_t v0_i = face.vertex_ids[e];

                const auto& v0 = grid_view->vertices[v0_i];
                chi_mesh::Vector3 n  = cell.faces[f].normal;

                chi_mesh::Vector3 v0p = interpolation_points[p] - v0;
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
              interpolation_points_ass_cell[p] = cell_local_index;
              interpolation_points_has_ass_cell[p] = true;
              chi_log.Log(LOG_ALLVERBOSE_2)
                << "Cell inter section found  " << p;
            }

          }//for each point

        }//if polyhedron
      }//for local cell

      //================================================== Upload node indices that
      //                                                   need mapping
      size_t num_interp_points = interpolation_points_ass_cell.size();
      for (size_t c=0; c<num_interp_points; c++)
      {
        if (not interpolation_points_has_ass_cell[c]) continue;

        uint64_t cell_local_index = interpolation_points_ass_cell[c];
        const auto& cell = grid_view->local_cells[cell_local_index];

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
        if (cell.Type() == chi_mesh::CellType::SLAB)
        {
          for (int i=0; i<2; i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(i);
            cfem_local_cells_needed_unmapped.push_back(cell_local_index);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_local_index);
          }

        }//if poly

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
        if (cell.Type() == chi_mesh::CellType::POLYGON)
        {
          size_t num_verts = cell.vertex_ids.size();
          for (size_t i=0; i<num_verts; i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(i);
            cfem_local_cells_needed_unmapped.push_back(cell_local_index);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_local_index);
          }

        }//if poly

          //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
        else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
        {
          size_t num_verts = cell.vertex_ids.size();
          for (size_t i=0; i<num_verts; i++)
          {
            cfem_local_nodes_needed_unmapped.push_back(i);
            cfem_local_cells_needed_unmapped.push_back(cell_local_index);
            pwld_local_nodes_needed_unmapped.push_back(i);
            pwld_local_cells_needed_unmapped.push_back(cell_local_index);
          }

        }//if polyh
      }//for associated cells
    }//if unique grid


    //Copies the latest developed references to the specific
    //field function context
    ff_context->cfem_local_nodes_needed_unmapped  = cfem_local_nodes_needed_unmapped;
    ff_context->cfem_local_cells_needed_unmapped  = cfem_local_cells_needed_unmapped;
    ff_context->pwld_local_nodes_needed_unmapped  = pwld_local_nodes_needed_unmapped;
    ff_context->pwld_local_cells_needed_unmapped  = pwld_local_cells_needed_unmapped;
    ff_context->interpolation_points_ass_cell     = interpolation_points_ass_cell;
    ff_context->interpolation_points_has_ass_cell = interpolation_points_has_ass_cell;
  }//for ff



  chi_log.Log(LOG_0VERBOSE_1) << "Finished initializing interpolator.";
}