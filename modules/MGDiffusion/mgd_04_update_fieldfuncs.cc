#include "mg_diffusion_solver.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

//###################################################################
/**Updates the field functions with the lates data.*/
void mg_diffusion::Solver::UpdateFieldFunctions()
{
  const auto& OneDOFPerNode = sdm_ptr->UNITARY_UNKNOWN_MANAGER;
  for (int g=0; g<num_groups; ++g)
  {
    std::vector<double> data_vector;
    sdm_ptr->LocalizePETScVector(x[g], data_vector, OneDOFPerNode);

    auto& ff = field_functions.at(g);
    ff->UpdateFieldVector(data_vector);
  }//for g
}