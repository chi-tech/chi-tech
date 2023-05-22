#ifndef CHITECH_FIELD_OPERATIONS_MULTIFIELD_H
#define CHITECH_FIELD_OPERATIONS_MULTIFIELD_H

#include "field_operation.h"
#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"
#include "ChiMath/Functions/function_dimA_to_dimB.h"

namespace chi_physics::field_operations
{

/**A field operation to manipulate a single field on the hand
* of a number of other fields.*/
class MultiFieldOperation : public FieldOperation
{
private:
  const size_t result_field_handle_;
  const std::vector<size_t> dependent_field_handles_;
  const size_t function_handle_;

  std::vector<unsigned int> dependent_field_ref_component_;
  std::vector<unsigned int> result_component_references_;

  std::shared_ptr<FieldFunctionGridBased> primary_ff_;
  std::vector<std::shared_ptr<const FieldFunction>> dependent_ffs_;

  std::shared_ptr<const chi_math::FunctionDimAToDimB> function_ptr_;
public:
  static chi_objects::InputParameters GetInputParameters();

  explicit MultiFieldOperation(const chi_objects::InputParameters& params);

  void Execute() override;
};

}

#endif // CHITECH_FIELD_OPERATIONS_MULTIFIELD_H
