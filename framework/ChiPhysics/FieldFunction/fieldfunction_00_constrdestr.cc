#include <utility>

#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_physics
{

FieldFunction::FieldFunction(std::string text_name,
                             chi_math::SMDPtr &sdm_ptr,
                             chi_math::Unknown unknown) :
  text_name_(std::move(text_name)),
  sdm_(sdm_ptr),
  unknown_(std::move(unknown)),
  unknown_manager_({unknown})
{
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(unknown_manager_);
  field_vector_.assign(num_local_dofs, 0.0);
}

FieldFunction::FieldFunction(std::string text_name,
                             chi_math::SMDPtr &sdm_ptr,
                             chi_math::Unknown unknown,
                             std::vector<double>  field_vector) :
  text_name_(std::move(text_name)),
  sdm_(sdm_ptr),
  unknown_(std::move(unknown)),
  field_vector_(std::move(field_vector)),
  unknown_manager_({unknown})
{
  const std::string fname = __FUNCTION__;
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(unknown_manager_);
  if (field_vector.size() != num_local_dofs)
    throw std::logic_error(fname + ": Constructor initialized with incompatible "
                                   "size field vector.");

  field_vector_ = field_vector;
}

FieldFunction::FieldFunction(std::string text_name,
                             chi_math::SMDPtr &sdm_ptr,
                             chi_math::Unknown unknown,
                             double field_value) :
  text_name_(std::move(text_name)),
  sdm_(sdm_ptr),
  unknown_(std::move(unknown)),
  unknown_manager_({unknown})
{
  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(unknown_manager_);
  field_vector_.assign(num_local_dofs, field_value);
}

FieldFunction::FieldFunction(std::string text_name,
                             chi_math::SMDPtr &sdm_ptr,
                             chi_math::Unknown unknown,
                             const std::vector<double>&  field_component_value) :
  text_name_(std::move(text_name)),
  sdm_(sdm_ptr),
  unknown_(std::move(unknown)),
  unknown_manager_({unknown})
{
  const std::string fname = __FUNCTION__;
  if (field_component_value.size() != unknown.num_components_)
    throw std::logic_error(fname + ": Constructor initialized with incompatible "
                                   "number of component values.");

  const size_t num_local_dofs = sdm_->GetNumLocalDOFs(unknown_manager_);
  field_vector_.assign(num_local_dofs, 0.0);

  const auto& grid = sdm_ptr->ref_grid_;
  for (const auto& cell : grid.local_cells)
  {
    const size_t num_nodes = sdm_->GetCellNumNodes(cell);
    for (size_t n=0; n<num_nodes; ++n)
      for (size_t c=0; c < unknown_.num_components_; ++c)
      {
        const auto nmap = sdm_->MapDOFLocal(cell, n, unknown_manager_, 0, c);

        if (nmap >= 0) field_vector_[nmap] = field_component_value[c];
      }
  }//for cell
}


}//namespace chi_physics