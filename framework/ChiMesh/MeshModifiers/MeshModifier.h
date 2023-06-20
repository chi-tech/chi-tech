#ifndef CHITECH_MESHMODIFIER_H
#define CHITECH_MESHMODIFIER_H

#include "ChiObject/chi_object.h"

namespace chi_mesh
{

/**Base class for mesh modifiers*/
class MeshModifier : public ChiObject
{
public:
  explicit MeshModifier(const chi_objects::InputParameters& params);

  virtual void Apply() = 0;
  virtual ~MeshModifier() = default;
protected:
  MeshModifier() = default;
};

}

#endif // CHITECH_MESHMODIFIER_H
