#ifndef CHI_MATH_UNKNOWN_MANAGER_H
#define CHI_MATH_UNKNOWN_MANAGER_H

#include <vector>
#include <string>
#include <stdexcept>

namespace chi_math
{

//#################################################################
/**Different types of variables.*/
enum class NodalVariableType
{
  SCALAR   = 1,
  VECTOR_2 = 2,
  VECTOR_3 = 3,
  VECTOR_N = 4,
  TENSOR   = 5
};

/**Nodal variable storage format.*/
enum class NodalStorageType
{
  NODAL = 1,
  BLOCK = 2
};

//###################################################################
/**General object for the management of unknowns in mesh-based
 * mathematical model.*/
class NodalVariableStructure
{
private:
  //#################################################################
  /**Basic class for an variable.*/
  class Variable
  {
  public:
    const NodalVariableType type;
    const unsigned int num_components;
    const unsigned int map_begin;
    std::string text_name;
    std::vector<std::string> component_text_names;
    std::vector<int> num_off_block_connections;

  public:
    explicit Variable(NodalVariableType in_type,
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
        case NodalVariableType::SCALAR:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=1"
                                    " for a SCALAR unknown.");
          map_value = 0;
          break;
        case NodalVariableType::VECTOR_2:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=2"
                                    " for a VECTOR_2 unknown.");
          map_value = map_begin + component_number;
          break;
        case NodalVariableType::VECTOR_3:
          if (component_number >= num_components)
            throw std::out_of_range("Attempting to access component >=3"
                                    " for a VECTOR_3 unknown.");
          map_value = map_begin + component_number;
          break;
        case NodalVariableType::VECTOR_N:
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
  std::vector<Variable> unknowns;
  const NodalStorageType dof_storage_type;

  explicit NodalVariableStructure(NodalStorageType in_storage_type=
                                  NodalStorageType::NODAL) :
    dof_storage_type(in_storage_type)
  {}

  unsigned int AddVariable(NodalVariableType unk_type, unsigned int dimension= 0);

  unsigned int MapVariable(unsigned int unknown_id, unsigned int component = 0);

  unsigned int GetTotalVariableStructureSize();

  void SetVariableNumOffBlockConnections(unsigned int unknown_id,
                                         int num_conn);
  void SetVariableComponentNumOffBlockConnections(unsigned int unknown_id,
                                                  unsigned int component,
                                                  int num_conn);
  void SetVariableTextName(unsigned int unknown_id,
                           const std::string& in_text_name);
  void SetVariableComponentTextName(unsigned int unknown_id,
                                    unsigned int component,
                                    const std::string& in_text_name);

  ~NodalVariableStructure() = default;
};







}

#endif