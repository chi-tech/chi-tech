#ifndef _chi_math_unknown_manager_h
#define _chi_math_unknown_manager_h

#include <vector>
#include <string>

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

  enum class DOFStorageType
  {
    NODAL = 1,
    BLOCK = 2
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
    const unsigned int num_components;
    std::string text_name;
    std::vector<std::string> component_text_names;
    std::vector<int> num_off_block_connections;

  public:
    explicit Unknown(UnknownType in_type, int num_comps) :
      type(in_type),
      num_components(num_comps)
    {
      component_text_names.resize(num_comps);
      num_off_block_connections.resize(num_comps,0);
    }

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
      Unknown(UnknownType::SCALAR,1),
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
  public:
    Vector2Unknown(unsigned int in_map_begin) :
      Unknown(UnknownType::VECTOR_2,2),
      map_begin(in_map_begin)
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
  public:
    Vector3Unknown(unsigned int in_map_begin) :
      Unknown(UnknownType::VECTOR_3,3),
      map_begin(in_map_begin)
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
  public:
    VectorNUnknown(unsigned int in_map_begin, unsigned int dimension) :
      Unknown(UnknownType::VECTOR_N,dimension),
      map_begin(in_map_begin)
    {}

    unsigned int GetMap(unsigned int component_0) override
    {return map_begin + component_0;}

    unsigned int GetMapEnd() override
    {return map_begin + num_components - 1;}
  };


public:
  std::vector<Unknown*> unknowns;
  const DOFStorageType dof_storage_type;

  explicit UnknownManager(DOFStorageType in_storage_type=DOFStorageType::NODAL) :
    dof_storage_type(in_storage_type)
  {}

  unsigned int AddUnknown(UnknownType unk_type, unsigned int dimension=0);

  unsigned int MapUnknown(unsigned int unknown_id, unsigned int component = 0);

  unsigned int GetTotalUnknownSize();

  void SetUnknownNumOffBlockConnections(unsigned int unknown_id,
                                        int num_conn);
  void SetComponentNumOffBlockConnections(unsigned int unknown_id,
                                          unsigned int component,
                                          int num_conn);
  void SetUnknownTextName(unsigned int unknown_id,
                          const std::string& in_text_name);
  void SetUnknownComponentTextName(unsigned int unknown_id,
                                   unsigned int component,
                                   const std::string& in_text_name);

  ~UnknownManager()
  {
    for (auto unknown : unknowns)
      delete unknown;
  }
};







}

#endif