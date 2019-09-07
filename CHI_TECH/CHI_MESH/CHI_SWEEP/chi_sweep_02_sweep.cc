#include "../chi_mesh.h"
#include "chi_sweep.h"

#include "../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../CHI_CELL/cell_polyhedron.h"
#include "../CHI_CELL/cell_polygon.h"

#include "../CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../CHI_VOLUMEMESHER/chi_volumemesher.h"

#include "../../CHI_GRAPH/chi_graph.h"

#include "chi_SPDS.h"
#include "chi_angleaggregation.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include "../../ChiTimer/chi_timer.h"

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

extern double chi_global_timings[20];

//###################################################################
/**Performs a sweep for the given angle aggregation.*/
//void chi_mesh::SweepManagement::
//Sweep(AngleAggregation* angle_agg,
//      chi_mesh::SweepManagement::SweepChunk* sweep_chunk)
//{
//  ChiTimer t16_sweeptime; t16_sweeptime.Reset();
//  //================================================== Loop over AngleSetGroups
//  // For 3D geometry this will be 8, one for each octant.
//  // For 2D geometry this will be 4, one for each quadrant.
//  // For 1D geometry this will be 2, one for left and one for right
//  bool completion_status = FLAG_FINISHED;
//  bool first_group = true;
//  while ((completion_status == FLAG_NOT_FINISHED) || (first_group))
//  {
//    completion_status = FLAG_FINISHED;
//    first_group=false;
//    for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
//    {
//      bool status = angle_agg->angle_set_groups[q]->
//        AngleSetGroupAdvance(sweep_chunk, q);
//      completion_status = completion_status && status;
//    }
//  }
//
//  //================================================== Reset all
//  for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
//  {
//    angle_agg->angle_set_groups[q]->ResetSweep();
//  }
//  chi_global_timings[16] += t16_sweeptime.GetTime()/1000.0;
//  chi_global_timings[17] += 1.0;
//
//
//}


