#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"

//#########################################################
/**Makes a ghosted version of the field vector.*/
std::vector<double> chi_physics::FieldFunction::GetGhostedFieldVector() const
{
  const size_t num_local_dofs = m_sdm->GetNumLocalDOFs(m_unknown_manager);
  const size_t num_globl_dofs = m_sdm->GetNumGlobalDOFs(m_unknown_manager);
  const std::vector<int64_t> ghost_ids =
    m_sdm->GetGhostDOFIndices(m_unknown_manager);

  chi_math::VectorGhostCommunicator vgc(num_local_dofs,
                                        num_globl_dofs,
                                        ghost_ids,
                                        MPI_COMM_WORLD);
  std::vector<double> field_wg = vgc.MakeGhostedVector(FieldVectorRead());

  vgc.CommunicateGhostEntries(field_wg);

  return field_wg;
}