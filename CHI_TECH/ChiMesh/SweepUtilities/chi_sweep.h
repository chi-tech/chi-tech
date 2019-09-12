#ifndef _chi_sweep_h
#define _chi_sweep_h

#include "../chi_mesh.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
namespace chi_mesh
{
namespace SweepManagement
{
  struct STDG;   ///< Global Sweep Plane Ordering
  struct SPLS;   ///< Sweep Plane Local Subgrid
  class  FLUDS;  ///< Flux Data Structure
  struct SPDS;   ///< Sweep Plane Data Structure

  class  SweepBuffer;
  class AngleSet;
  class AngleSetGroup;
  class  AngleAggregation;

  class SweepChunk;

  class SweepScheduler;

  //01
  SPDS* CreateSweepOrder(double polar, double azimuthal,
                         chi_mesh::MeshContinuum *vol_continuum,
                         int number_of_groups);

  //02
//  void Sweep(SPDS* spds,
//             SweepChunk* sweep_chunk=NULL,
//             bool communicate=true);
//  void Sweep(AngleAggregation* angle_agg,
//             SweepChunk* sweep_chunk=NULL);
  //03
  void PrintSweepOrdering(SPDS* sweep_order,
                          MeshContinuum* vol_continuum);
}
}

#include "chi_sweepchunk_base.h"

#endif
