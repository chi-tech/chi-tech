#include "fieldfunction_gridbased.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

// #########################################################
/**Makes a ghosted version of the field vector.*/
std::vector<double> chi_physics::FieldFunctionGridBased::GetGhostedFieldVector() const
{
  // const size_t num_local_dofs = sdm_->GetNumLocalDOFs(unknown_manager_);
  // const size_t num_globl_dofs = sdm_->GetNumGlobalDOFs(unknown_manager_);
  // const std::vector<int64_t> ghost_ids =
  //   sdm_->GetGhostDOFIndices(unknown_manager_);
  //
  // chi_math::VectorGhostCommunicator vgc(
  //   num_local_dofs, num_globl_dofs, ghost_ids, MPI_COMM_WORLD);
  // std::vector<double> field_wg = vgc.MakeGhostedVector(FieldVectorRead());

  // vgc.CommunicateGhostEntries(field_wg);

  const auto& vgc = *vector_ghost_communicator_;

   std::vector<double> field_wg = vgc.MakeGhostedVector(FieldVectorRead());

   vgc.CommunicateGhostEntries(field_wg);

  return field_wg;
}