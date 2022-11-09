#ifndef CHITECH_FIELDFUNCTION2_H
#define CHITECH_FIELDFUNCTION2_H

#include <string>
#include <memory>
#include <vector>

#include "ChiMath/UnknownManager/unknown_manager.h"

#include <vtkUnstructuredGrid.h>

namespace chi_mesh
{
  class Cell;
}

//######################################################### Forward declarations
namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<const MeshContinuum> MeshContinuumConstPtr;
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
  std::string                     m_text_name;
  chi_math::SMDPtr                m_sdm;
  chi_math::Unknown               m_unknown;

  std::vector<double>             m_field_vector;

  chi_math::UnknownManager        m_unknown_manager;

public:
  FieldFunction2(std::string                      text_name,
                 chi_math::SMDPtr&                sdm_ptr,
                 chi_math::Unknown                unknown);

  FieldFunction2(std::string                      text_name,
                 chi_math::SMDPtr&                sdm_ptr,
                 chi_math::Unknown                unknown,
                 std::vector<double>              field_vector);

  FieldFunction2(std::string                      text_name,
                 chi_math::SMDPtr&                sdm_ptr,
                 chi_math::Unknown                unknown,
                 double                           field_value);

  FieldFunction2(std::string                      text_name,
                 chi_math::SMDPtr&                sdm_ptr,
                 chi_math::Unknown                unknown,
                 const std::vector<double>&       field_component_value);

  //01 Updates
  void UpdateFieldVector(const std::vector<double>& field_vector);
  const std::vector<double>& FieldVector() const;

  //03 Export VTK
  void ExportToVTK(const std::string& file_base_name) const;

private:
  static void UploadCellGeometry(const chi_mesh::MeshContinuum& grid,
                          const chi_mesh::Cell &cell,
                          int64_t& node_counter,
                          vtkNew<vtkPoints>& points,
                          vtkNew<vtkUnstructuredGrid> &ugrid);

public:
  std::vector<double> GetGhostedFieldVector() const;
};

}//namespace chi_physics

#endif //CHITECH_FIELDFUNCTION2_H
