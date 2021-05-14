#include "lbs_linear_boltzmann_solver.h"
#include "ChiConsole/chi_console.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include <iomanip>

extern ChiConsole&  chi_console;
extern ChiLog&     chi_log;
extern ChiMPI&      chi_mpi;

#include "ChiTimer/chi_timer.h"

extern ChiTimer chi_program_timer;

//###################################################################
/**Initializes fluds data structures.*/
void LinearBoltzmann::Solver::InitFluxDataStructures(LBSGroupset& groupset)
{
  //================================================== Initialize grid vertex
  //                                                   mapping
  // This is crucial information for the transfer of
  // information during a sweep. And is used my all the
  // flux data structures.
  grid_nodal_mappings.clear();
  grid_nodal_mappings.reserve(grid->local_cells.size());
  for (auto& cell : grid->local_cells)
  {
    using namespace chi_mesh::sweep_management;
    CellFaceNodalMapping cell_nodal_mapping;

    cell_nodal_mapping.reserve(cell.faces.size());
    for (auto& face : cell.faces)
    {
      std::vector<short> face_nodal_mapping;
      int ass_face = -1;
      if (face.has_neighbor and face.IsNeighborLocal(*grid))
      {
        grid->FindAssociatedVertices(face,face_nodal_mapping);
        ass_face = face.GetNeighborAssociatedFace(*grid);
      }

      cell_nodal_mapping.emplace_back(ass_face,face_nodal_mapping);
    }//for f

    grid_nodal_mappings.push_back(std::move(cell_nodal_mapping));
  }//for local cell

  //================================================== Angle Aggregation
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher& mesher = *handler->volume_mesher;

  if ( options.geometry_type == GeometryType::ONED_SLAB or
       options.geometry_type == GeometryType::TWOD_CARTESIAN or
       (typeid(mesher) == typeid(chi_mesh::VolumeMesherExtruder)))
  {
    switch (groupset.angleagg_method)
    {
      case AngleAggregationType::SINGLE:
        InitAngleAggSingle(groupset); break;
      case AngleAggregationType::POLAR:
        InitAngleAggPolar(groupset); break;
      default:
        throw std::logic_error(std::string(__FUNCTION__) +
                               " Invalid angle aggregation type.");
    }//switch on method
  }//if aggregatable
  else if (options.geometry_type == GeometryType::ONED_SPHERICAL ||
           options.geometry_type == GeometryType::TWOD_CYLINDRICAL)
  {
    switch (groupset.angleagg_method)
    {
      case AngleAggregationType::SINGLE:
        InitAngleAggSingle(groupset); break;
      case AngleAggregationType::AZIMUTHAL:
        InitAngleAggAzimuthal(groupset); break;
      default:
        throw std::logic_error(std::string(__FUNCTION__) +
                               " Invalid angle aggregation type.");
    }//switch on method
  }//if aggregatable
  else
    InitAngleAggSingle(groupset);

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Initialized Angle Aggregation.   "
    << "         Process memory = "
    << std::setprecision(3) << chi_console.GetMemoryUsageInMB()
    << " MB.";


  MPI_Barrier(MPI_COMM_WORLD);
}
