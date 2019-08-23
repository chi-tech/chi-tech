#include "chi_angleset.h"
#include "chi_SPDS.h"
#include "chi_sweepbuffer.h"

#include <chi_mpi.h>

extern CHI_MPI chi_mpi;

#include "../../CHI_TIMER/chi_timer.h"


//###################################################################
/**This function advances the work stages of an angleset.*/
bool chi_mesh::SweepManagement::AngleSet::
AngleSetAdvance(chi_mesh::SweepManagement::SweepChunk *sweep_chunk,
                int angle_set_num)
{
  //================================================== Prevent reexecution
  if (executed)
    return FLAG_FINISHED;

  //================================================== Check sweep buffer
  sweep_buffer.CheckInitialized();

  //================================================== Check all predecessor
  //                                                   locations sent data
  if (!sweep_buffer.CheckUpstreamPsiAvailable(angle_set_num))
    return FLAG_NOT_FINISHED;

  //================================================== Resize the angle set's
  //                                                   FLUDS at workstage start
  sweep_buffer.InitializeBuffers();

  //================================================== Execute chunk
  sweep_chunk->Sweep(this);

  //================================================== Send outgoing psi
  sweep_buffer.SendDownstreamPsi(angle_set_num);

  //================================================== Release memory
  sweep_buffer.ClearReceiveBuffers();

  executed = true;
  return FLAG_FINISHED;
}