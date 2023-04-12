#ifndef CHITECH_OBJECT_MAKER_H
#define CHITECH_OBJECT_MAKER_H

#include <memory>
#include "ChiParameters/input_parameters.h"
#include "ChiParameters/parameter_block.h"
#include "chi_object.h"

/**Small utility macro for joining two words.*/
#define ChiObjectMakerJoinWordsA(x, y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define ChiObjectJoinWordsB(x, y) ChiObjectMakerJoinWordsA(x, y)

/**Macro for registering an object within the ObjectMaker singleton.
 * Example:
 *
 * \note Remember to include the header "ChiObject/object_maker.h"*/
#define RegisterChiObject(namespace_name, object_name)                         \
  static char ChiObjectJoinWordsB(unique_var_name_object_##object_name##_,     \
                                  __COUNTER__) =                               \
    chi_objects::ObjectMaker::AddObjectToRegistry<object_name,                 \
                                                  chi_objects::ChiObject>(     \
      #namespace_name, #object_name)

namespace chi_objects
{

class ChiObject;

class ObjectMaker
{
private:
  using ObjectPtr = std::shared_ptr<ChiObject>;

  using ObjectGetInParamsFunc = InputParameters (*)();
  using ObjectConstructorFunc = ObjectPtr (*)(const InputParameters&);

  struct ObjectRegistryEntry
  {
    ObjectGetInParamsFunc get_in_params_func = nullptr;
    ObjectConstructorFunc constructor_func = nullptr;
  };

  std::map<std::string, ObjectRegistryEntry> object_registry_;

private:
  ObjectMaker() = default;

public:
  ObjectMaker(const ObjectMaker&) = delete;
  ObjectMaker(const ObjectMaker&&) = delete;
  ObjectMaker& operator=(const ObjectMaker&) = delete;

  static ObjectMaker& GetInstance() noexcept;

  const std::map<std::string, ObjectRegistryEntry>& Registry() const;

private:
  template <typename T>
  static chi_objects::InputParameters CallGetInputParamsFunction()
  {
    return T::GetInputParameters();
  }

  template <typename T, typename base_T>
  static std::shared_ptr<base_T>
  CallObjectConstructor(const InputParameters& params)
  {
    return std::make_shared<T>(params);
  }

public:
  template <typename T, typename base_T>
  static char AddObjectToRegistry(const std::string& namespace_name,
                                  const std::string& object_name)
  {
    const std::string name = namespace_name + "::" + object_name;

    auto& object_maker = GetInstance();

    // Check if the function name is already there
    if (object_maker.object_registry_.count(name) > 0)
    {
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Attempted "
                             "to register Object \"" +
                             name +
                             "\" but an object with the same name is"
                             " already registered.");
    }

    ObjectRegistryEntry reg_entry;
    reg_entry.get_in_params_func = &CallGetInputParamsFunction<T>;
    reg_entry.constructor_func = &CallObjectConstructor<T, base_T>;
    object_maker.object_registry_.insert(std::make_pair(name, reg_entry));

    return 0;
  }

  size_t MakeObject(const ParameterBlock& params) const;
  size_t MakeObjectType(const std::string& type,
                        const ParameterBlock& params) const;

  /**Dumps the object registry to stdout.*/
  void DumpRegister() const;
};

} // namespace chi_objects

#endif // CHITECH_OBJECT_MAKER_H
