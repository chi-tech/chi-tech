#ifndef chi_runtime_h
#define chi_runtime_h

#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

namespace chi_mesh
{
  class MeshHandler;
}


/**General utilities in ChiTech*/
class chi
{
public:
  static std::vector<std::shared_ptr<chi_mesh::MeshHandler>>  meshhandler_stack;
  static int         current_mesh_handler;

  class run_time
  {
  public:
    static bool        termination_posted;
    static std::string input_file_name;
    static bool        sim_option_interactive;
    static bool        allow_petsc_error_handler;

  private:
    static void ParseArguments(int argc, char** argv);

  public:
    static int  RunInteractive(int argc, char** argv);
    static int  RunBatch(int argc, char** argv);
    static int  Initialize(int argc, char** argv);
    static void Finalize();
  };


public:
  /**Attempts to retrieve an object of base-type `shared_ptr<T>` at the given
   * handle. It then attempts to cast it to type `shared_ptr<R>` and, if
   * successful, will return a reference of type R&.*/
  template<class R,class T>
  static R& GetStackItem(std::vector<std::shared_ptr<T>>& stack,
                        const size_t handle,
                        const std::string& calling_function_name="Unknown")
  {
    try
    {
      std::shared_ptr<T>& item = stack.at(handle);
      std::shared_ptr<R> ret_item = std::dynamic_pointer_cast<R>(item);
      if (not ret_item)
        throw std::logic_error("chi::GetStackItem: Invalid return type used. "
                               "Calling function: " + calling_function_name);
      return *ret_item;
    }
    catch (const std::out_of_range& oor)
    {
      throw std::logic_error("chi::GetStackItem: Invalid handle used. "
                             "Calling function: " + calling_function_name);
    }
  }
  chi() = delete;
};

#endif