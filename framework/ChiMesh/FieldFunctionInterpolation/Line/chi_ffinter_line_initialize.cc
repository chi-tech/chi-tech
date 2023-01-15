#include "chi_ffinter_line.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Initializes the data structures necessary for interpolation. This is
 * independent of the physics and hence is a routine on its own.*/
void chi_mesh::FieldFunctionInterpolationLine::
  Initialize()
{
  chi::log.Log0Verbose1() << "Initializing line interpolator.";
  //================================================== Check for empty FF-list
  if (field_functions.empty())
    throw std::logic_error("Unassigned field function in line "
                           "field function interpolator.");

  //================================================== Create points;
  const chi_mesh::Vector3 vif = pf - pi;
  delta_d = vif.Norm()/(number_of_points-1);

  const auto omega = vif.Normalized();

  interpolation_points.push_back(pi);
  for (int k=1; k<(number_of_points); k++)
    interpolation_points.push_back(pi + omega*delta_d*k);

  //====================================================== Loop over contexts
  const size_t num_ff = field_functions.size();
  for (size_t ff=0; ff<num_ff; ff++)
  {
    ff_contexts.emplace_back();
    auto& ff_context = ff_contexts.back();

    ff_context.ref_ff = field_functions[ff];
    const auto& sdm = ff_context.ref_ff->SDM();
    const auto& grid = *sdm.ref_grid;

    ff_context.interpolation_points_ass_cell.assign(number_of_points,0);
    ff_context.interpolation_points_has_ass_cell.assign(number_of_points,false);

    //================================================== Find a home for each
    //                                                   point
    for (const auto& cell : grid.local_cells)
    {
      for (int p=0; p<number_of_points; p++)
      {
        const auto& point = interpolation_points[p];
        if (grid.CheckPointInsideCell(cell, point))
        {
          ff_context.interpolation_points_ass_cell[p] = cell.local_id;
          ff_context.interpolation_points_has_ass_cell[p] = true;
        }
      }//for point p
    }//for cell
  }//for ff

  chi::log.Log0Verbose1() << "Finished initializing interpolator.";
}