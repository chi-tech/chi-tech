#ifndef CHITECH_FIELDFUNCTION_GRIDBASED_H
#define CHITECH_FIELDFUNCTION_GRIDBASED_H

#include "fieldfunction.h"
#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "ChiMesh/chi_mesh.h"

#include "vtkUnstructuredGrid.h"
#include <petscksp.h>

#include <string>
#include <memory>
#include <vector>
#include <utility>

// ######################################################### Forward decls
namespace chi_math
{
class SpatialDiscretization;
typedef std::shared_ptr<SpatialDiscretization> SMDPtr;
} // namespace chi_math

namespace chi_physics
{

// ################################################################### Class def
/***/
class FieldFunctionGridBased : public FieldFunction
{
public:
  typedef std::pair<chi_mesh::Vector3, chi_mesh::Vector3> BoundingBox;
  typedef std::shared_ptr<chi_math::VectorGhostCommunicator> VectorGhostCommPtr;

protected:
  chi_math::SMDPtr sdm_;
  std::vector<double> field_vector_;

private:
  const BoundingBox local_grid_bounding_box_;
  VectorGhostCommPtr vector_ghost_communicator_ = nullptr;

public:
  /**Returns required input parameters.*/
  static chi::InputParameters GetInputParameters();

  /**ObjectMaker based constructor.*/
  explicit FieldFunctionGridBased(const chi::InputParameters& params);

  /**Creates a field function, filling it with zeros.*/
  FieldFunctionGridBased(const std::string& text_name,
                         chi_math::SMDPtr& sdm_ptr,
                         chi_math::Unknown unknown);

  /**Creates a field function with an associated field vector.
   * The field's data vector is set to the incoming field vector.*/
  FieldFunctionGridBased(const std::string& text_name,
                         chi_math::SMDPtr& sdm_ptr,
                         chi_math::Unknown unknown,
                         std::vector<double> field_vector);

  /**Creates a field function where all the values are assigned to
   * the single supplied value.*/
  FieldFunctionGridBased(const std::string& text_name,
                         chi_math::SMDPtr& sdm_ptr,
                         chi_math::Unknown unknown,
                         double field_value);

  virtual ~FieldFunctionGridBased() = default;

private:
  VectorGhostCommPtr MakeGhostCommunicator();

public:
  // Getters
  const chi_math::SpatialDiscretization& SDM() const { return *sdm_; }
  const std::vector<double>& FieldVectorRead() const { return field_vector_; }
  std::vector<double>& FieldVector() { return field_vector_; }

  // 01 Updates
  void UpdateFieldVector(const std::vector<double>& field_vector);
  void UpdateFieldVector(const Vec& field_vector);

  // 03 Export VTK
  typedef std::vector<std::shared_ptr<const FieldFunctionGridBased>> FFList;
  static void ExportMultipleToVTK(const std::string& file_base_name,
                                  const FFList& ff_list);

public:
  // 04 Utils
  std::vector<double> GetGhostedFieldVector() const;

  // 05 Point Values
  /**\brief Returns the component values at requested point.*/
  virtual std::vector<double>
  GetPointValue(const chi_mesh::Vector3& point) const;

  double Evaluate(const chi_mesh::Cell& cell,
                  const chi_mesh::Vector3& position,
                  unsigned int component) const override;
};

} // namespace chi_physics

#endif // CHITECH_FIELDFUNCTION_GRIDBASED_H
