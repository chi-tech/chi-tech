#ifndef CHI_MATH_UNKNOWN_MANAGER_H
#define CHI_MATH_UNKNOWN_MANAGER_H

#include <vector>
#include <string>
#include <stdexcept>

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

/**Degree-of-freedom storage format.*/
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
  /**Basic class for an unknown.*/
  class Unknown
  {
  public:
    const UnknownType type;
    const unsigned int num_components;
    const unsigned int map_begin;
    std::string text_name;
    std::vector<std::string> component_text_names;
    std::vector<int> num_off_block_connections;

  public:
    explicit Unknown(UnknownType in_type,
                     unsigned int in_num_components=1,
                     unsigned int in_map_begin=0) :
      type(in_type),
      num_components(in_num_components),
      map_begin(in_map_begin)
    {
      component_text_names.resize(in_num_components,std::string());
      num_off_block_connections.resize(in_num_components, 0);
    }

    unsigned int GetMap(unsigned int component_number=0) const
    {
      unsigned int map_value = 0;
      switch (type)
      {
        case UnknownType::SCALAR:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=1"
                                    " for a SCALAR unknown.");
          map_value = 0;
          break;
        case UnknownType::VECTOR_2:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=2"
                                    " for a VECTOR_2 unknown.");
          map_value = map_begin + component_number;
          break;
        case UnknownType::VECTOR_3:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=3"
                                    " for a VECTOR_3 unknown.");
          map_value = map_begin + component_number;
          break;
        case UnknownType::VECTOR_N:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >="+
                                    std::to_string(num_components)+
                                    " for a VECTOR_2 unknown.");
          map_value = map_begin + component_number;
          break;
        default:
          break;
      }

      return map_value;
    }
    unsigned int GetMapEnd() const {return map_begin + num_components - 1;}
  };


public:
  std::vector<Unknown> unknowns;
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

  ~UnknownManager() = default;
};







}

#endif