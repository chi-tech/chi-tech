#include <utility>

#include "fieldfunction2.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_physics
{

FieldFunction2::FieldFunction2(std::string text_name,
                               chi_mesh::MeshContinuumPtr &grid_ptr,
                               chi_math::SMDPtr &sdm_ptr,
                               chi_math::Unknown unknown) :
                               m_text_name(std::move(text_name)),
                               m_grid(grid_ptr),
                               m_sdm(sdm_ptr),
                               m_unknown(std::move(unknown)),
                               m_unknown_manager({unknown})
{
  const size_t num_local_dofs = m_sdm->GetNumLocalDOFs(m_unknown_manager);
  m_field_vector.assign(num_local_dofs, 0.0);
}

FieldFunction2::FieldFunction2(std::string text_name,
                               chi_mesh::MeshContinuumPtr &grid_ptr,
                               chi_math::SMDPtr &sdm_ptr,
                               chi_math::Unknown unknown,
                               std::vector<double>  field_vector) :
  m_text_name(std::move(text_name)),
  m_grid(grid_ptr),
  m_sdm(sdm_ptr),
  m_unknown(std::move(unknown)),
  m_field_vector(std::move(field_vector)),
  m_unknown_manager({unknown})
{
  const std::string fname = __FUNCTION__;
  const size_t num_local_dofs = m_sdm->GetNumLocalDOFs(m_unknown_manager);
  if (field_vector.size() != num_local_dofs)
    throw std::logic_error(fname + ": Constructor initialized with incompatible "
                                   "size field vector.");

  m_field_vector = field_vector;
}

FieldFunction2::FieldFunction2(std::string text_name,
                               chi_mesh::MeshContinuumPtr &grid_ptr,
                               chi_math::SMDPtr &sdm_ptr,
                               chi_math::Unknown unknown,
                               double field_value) :
  m_text_name(std::move(text_name)),
  m_grid(grid_ptr),
  m_sdm(sdm_ptr),
  m_unknown(std::move(unknown)),
  m_unknown_manager({unknown})
{
  const size_t num_local_dofs = m_sdm->GetNumLocalDOFs(m_unknown_manager);
  m_field_vector.assign(num_local_dofs, field_value);
}

FieldFunction2::FieldFunction2(std::string text_name,
                               chi_mesh::MeshContinuumPtr &grid_ptr,
                               chi_math::SMDPtr &sdm_ptr,
                               chi_math::Unknown unknown,
                               const std::vector<double>&  field_component_value) :
  m_text_name(std::move(text_name)),
  m_grid(grid_ptr),
  m_sdm(sdm_ptr),
  m_unknown(std::move(unknown)),
  m_unknown_manager({unknown})
{
  const std::string fname = __FUNCTION__;
  if (field_component_value.size() != unknown.num_components)
    throw std::logic_error(fname + ": Constructor initialized with incompatible "
                                   "number of component values.");

  const size_t num_local_dofs = m_sdm->GetNumLocalDOFs(m_unknown_manager);
  m_field_vector.assign(num_local_dofs, 0.0);

  for (const auto& cell : m_grid->local_cells)
  {
    const size_t num_nodes = m_sdm->GetCellNumNodes(cell);
    for (size_t n=0; n<num_nodes; ++n)
      for (size_t c=0; c<m_unknown.num_components; ++c)
      {
        const auto nmap = m_sdm->MapDOFLocal(cell, n, m_unknown_manager, 0, c);

        if (nmap >= 0) m_field_vector[nmap] = field_component_value[c];
      }
  }//for cell
}


}//namespace chi_physics