#include "chi_ffinter_volume.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiPhysics/FieldFunction2/fieldfunction2.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Initializes the volume field function interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Initialize()
{
  chi::log.Log0Verbose1() << "Initializing volume interpolator.";
  //================================================== Check grid available
  if (field_functions.empty())
    throw std::logic_error("Unassigned field function in volume field "
                           "function interpolator.");

  if (logical_volume == nullptr)
    throw std::logic_error("Unassigned logical volume in volume field function"
                           "interpolator.");

  const auto& grid = *field_functions.front()->SDM().ref_grid;

  //================================================== Find cells inside volume
  for (const auto& cell : grid.local_cells)
    if (logical_volume->Inside(cell.centroid))
      cell_local_ids_inside_logvol.push_back(cell.local_id);
}