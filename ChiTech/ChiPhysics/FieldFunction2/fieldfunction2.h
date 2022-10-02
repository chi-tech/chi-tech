#ifndef CHITECH_FIELDFUNCTION2_H
#define CHITECH_FIELDFUNCTION2_H

#include <string>
#include <memory>
#include <vector>

#include "ChiMath/UnknownManager/unknown_manager.h"

//######################################################### Forward declarations
namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<chi_mesh::MeshContinuum> MeshContinuumPtr;
}

namespace chi_math
{
  class SpatialDiscretization;
  typedef std::shared_ptr<SpatialDiscretization> SMDPtr;
}

namespace chi_physics
{

//################################################################### Class def
/***/
class FieldFunction2
{
protected:
  std::string                m_text_name;
  chi_mesh::MeshContinuumPtr m_grid;
  chi_math::SMDPtr           m_sdm;
  chi_math::Unknown          m_unknown;

  std::vector<double>        m_field_vector;

  chi_math::UnknownManager   m_unknown_manager;

public:
  FieldFunction2(std::string text_name,
                 chi_mesh::MeshContinuumPtr& grid_ptr,
                 chi_math::SMDPtr&           sdm_ptr,
                 chi_math::Unknown     unknown);
  FieldFunction2(std::string text_name,
                 chi_mesh::MeshContinuumPtr& grid_ptr,
                 chi_math::SMDPtr&           sdm_ptr,
                 chi_math::Unknown     unknown,
                 std::vector<double>   field_vector);

  FieldFunction2(std::string text_name,
                 chi_mesh::MeshContinuumPtr& grid_ptr,
                 chi_math::SMDPtr&           sdm_ptr,
                 chi_math::Unknown           unknown,
                 double                      field_value);

  FieldFunction2(std::string text_name,
                 chi_mesh::MeshContinuumPtr& grid_ptr,
                 chi_math::SMDPtr&           sdm_ptr,
                 chi_math::Unknown unknown,
                 const std::vector<double>&  field_component_value);

  void UpdateFieldVector(const std::vector<double>& field_vector);
};

}//namespace chi_physics

#endif //CHITECH_FIELDFUNCTION2_H
