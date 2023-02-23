#ifndef CHITECH_FIELDFUNCTION_H
#define CHITECH_FIELDFUNCTION_H

#include <string>
#include <memory>
#include <vector>

#include "ChiMath/UnknownManager/unknown_manager.h"

#include "vtkUnstructuredGrid.h"
#include <petscksp.h>

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
class FieldFunction
{
protected:
  std::string                     m_text_name;
  chi_math::SMDPtr                m_sdm;
  chi_math::Unknown               m_unknown;

  std::vector<double>             m_field_vector;

  chi_math::UnknownManager        m_unknown_manager;

public:
  FieldFunction(std::string                      text_name,
                chi_math::SMDPtr&                sdm_ptr,
                chi_math::Unknown                unknown);

  FieldFunction(std::string                      text_name,
                chi_math::SMDPtr&                sdm_ptr,
                chi_math::Unknown                unknown,
                std::vector<double>              field_vector);

  FieldFunction(std::string                      text_name,
                chi_math::SMDPtr&                sdm_ptr,
                chi_math::Unknown                unknown,
                double                           field_value);

  FieldFunction(std::string                      text_name,
                chi_math::SMDPtr&                sdm_ptr,
                chi_math::Unknown                unknown,
                const std::vector<double>&       field_component_value);

  //Getters
  const std::string& TextName() const {return m_text_name;}
  const chi_math::SpatialDiscretization& SDM() const {return *m_sdm;}
  const chi_math::Unknown& Unknown() const {return m_unknown;}
  const std::vector<double>& FieldVectorRead() const {return m_field_vector;}
  const chi_math::UnknownManager& UnkManager() const {return m_unknown_manager;}

  //01 Updates
  void UpdateFieldVector(const std::vector<double>& field_vector);
  void UpdateFieldVector(const Vec& field_vector);

  //03 Export VTK
  typedef std::vector<std::shared_ptr<const FieldFunction>> FFList;
  static
  void ExportMultipleToVTK(const std::string& file_base_name,
                           const FFList& ff_list);

public:
  //04 Utils
  std::vector<double> GetGhostedFieldVector() const;
};

}//namespace chi_physics

#endif //CHITECH_FIELDFUNCTION_H
