#ifndef CHITECH_FIELD_COPY_H
#define CHITECH_FIELD_COPY_H

#include "field_operation.h"
#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"

namespace chi_physics::field_operations
{

/**Field operaiton that copies components of one field to the
* components of another.*/
class FieldCopyOperation : public FieldOperation
{
private:
  const size_t to_field_handle_;
  const size_t from_field_handle_;

  std::vector<size_t> to_components_;
  std::vector<size_t> from_components_;

  std::shared_ptr<FieldFunctionGridBased> to_ff_;
  std::shared_ptr<const FieldFunctionGridBased> from_ff_;

public:
  static chi_objects::InputParameters GetInputParameters();

  explicit FieldCopyOperation(const chi_objects::InputParameters& params);

  void Execute() override;
};

}

#endif // CHITECH_FIELD_COPY_H
