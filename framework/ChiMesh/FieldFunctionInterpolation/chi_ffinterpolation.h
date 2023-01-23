#ifndef CHI_FFINTERPOLATION_H
#define CHI_FFINTERPOLATION_H

#include <memory>
#include <vector>

namespace chi_physics
{
  class FieldFunction;
  typedef std::shared_ptr<FieldFunction> FieldFunction2CPtr;
}

namespace chi_mesh::ff_interpolation
{
  enum class Type : int
  {
    SLICE = 1,
    LINE  = 2,
    VOLUME = 3,
    POINT = 4
  };

  enum class Operation : int
  {
    OP_SUM     = 10,
    OP_AVG     = 11,
    OP_MAX     = 12,
    OP_SUM_LUA = 13,
    OP_AVG_LUA = 14,
    OP_MAX_LUA = 15,
  };

  enum class Property : int
  {
    PROBEPOINT     = 0,
    SLICEPOINT     = 1,
    SLICENORMAL    = 2,
    SLICETANGENT   = 3,
    SLICEBINORM    = 4,
    OPERATION      = 5,
    LOGICAL_VOLUME = 8,

    ADD_FIELD_FUNCTION = 9,
    SET_FIELD_FUNCTIONS = 10,

    FIRSTPOINT     = 11,
    SECONDPOINT    = 12,
    NUMBEROFPOINTS = 13,
    CUSTOM_ARRAY   = 14,
  };
}

namespace chi_mesh
{

//###################################################################
/** Base class for field-function interpolation objects.*/
class FieldFunctionInterpolation
{
protected:
  ff_interpolation::Type m_type;
  unsigned int m_ref_component = 0;
public:
  std::vector<chi_physics::FieldFunction2CPtr> field_functions;

public:
  explicit
  FieldFunctionInterpolation(ff_interpolation::Type type) :
    m_type(type) {}

  ff_interpolation::Type Type() const {return m_type;}

  virtual void Initialize(){};
  virtual void Execute(){};

  virtual std::string GetDefaultFileBaseName() const = 0;
  virtual void ExportPython(std::string base_name) = 0;
};

}//namespace chi_mesh




#endif