#ifndef _chi_math_unknown_manager_h
#define _chi_math_unknown_manager_h

#include <vector>

namespace chi_math
{

  //#################################################################
  /**Different types of unknowns.*/
  enum class UnknownType
  {
    SCALAR   = 1,
    VECTOR_2 = 2,
    VECTOR_3 = 3,
    VECTOR_N = 4,
    TENSOR   = 5
  };

//###################################################################
/**General object for the management of unknowns in mesh-based
 * mathematical model.*/
class UnknownManager
{
private:
  //#################################################################
  /**Base class for unknown.*/
  class Unknown
  {
  public:
    const UnknownType type;

  public:
    explicit Unknown(UnknownType in_type) : type(in_type)
    {}

    virtual ~Unknown() = default;

    virtual unsigned int GetMap(unsigned int component_0) {return 0;}
    virtual unsigned int GetMapEnd() {return 0;}
  };

  //#################################################################
  /**Scalar unknown.*/
  class ScalarUnknown : public Unknown
  {
  private:
    const unsigned int map;
  public:
    explicit ScalarUnknown(unsigned int map_begin) :
      Unknown(UnknownType::SCALAR),
      map(map_begin)
    { }

    unsigned int GetMap(unsigned int component_0) override
    {return map;}

    unsigned int GetMapEnd() override
    {return map;}
  };

  //#################################################################
  /**2-component vector.*/
  class Vector2Unknown : public Unknown
  {
  private:
    const unsigned int map_begin;
    const unsigned int num_components;
  public:
    Vector2Unknown(unsigned int in_map_begin) :
      Unknown(UnknownType::VECTOR_2),
      map_begin(in_map_begin),
      num_components(2)
    {}

    unsigned int GetMap(unsigned int component_0) override
    {return map_begin + component_0;}

    unsigned int GetMapEnd() override
    {return map_begin + num_components - 1;}
  };

  //#################################################################
  /**3-component vector.*/
  class Vector3Unknown : public Unknown
  {
  private:
    const unsigned int map_begin;
    const unsigned int num_components;
  public:
    Vector3Unknown(unsigned int in_map_begin) :
      Unknown(UnknownType::VECTOR_3),
      map_begin(in_map_begin),
      num_components(3)
    {}

    unsigned int GetMap(unsigned int component_0) override
    {return map_begin + component_0;}

    unsigned int GetMapEnd() override
    {return map_begin + num_components - 1;}
  };

  //#################################################################
  /**Multi-component vector other that 2 dimensional or 3 dimensional.*/
  class VectorNUnknown : public Unknown
  {
  private:
    const unsigned int map_begin;
    const unsigned int num_components;
  public:
    VectorNUnknown(unsigned int in_map_begin, unsigned int dimension) :
      Unknown(UnknownType::VECTOR_N),
      map_begin(in_map_begin),
      num_components(dimension)
    {}

    unsigned int GetMap(unsigned int component_0) override
    {return map_begin + component_0;}

    unsigned int GetMapEnd() override
    {return map_begin + num_components - 1;}
  };


public:
  std::vector<Unknown*> unknowns;

  unsigned int AddUnknown(UnknownType unk_type, unsigned int dimension=0);

  unsigned int MapUnknown(unsigned int unknown_id, unsigned int component = 0);

  unsigned int GetTotalUnknownSize();

  ~UnknownManager()
  {
    for (auto unknown : unknowns)
      delete unknown;
  }
};







}

#endif