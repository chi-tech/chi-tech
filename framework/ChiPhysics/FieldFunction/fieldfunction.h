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
  std::string                     text_name_;
  chi_math::SMDPtr                sdm_;
  chi_math::Unknown               unknown_;

  std::vector<double>             field_vector_;

  chi_math::UnknownManager        unknown_manager_;

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
  const std::string& TextName() const {return text_name_;}
  const chi_math::SpatialDiscretization& SDM() const {return *sdm_;}
  const chi_math::Unknown& Unknown() const {return unknown_;}
  const std::vector<double>& FieldVectorRead() const {return field_vector_;}
  const chi_math::UnknownManager& UnkManager() const {return unknown_manager_;}

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
