#ifndef chi_runtime_h
#define chi_runtime_h

#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

namespace chi_mesh
{
  class MeshHandler;
  typedef std::shared_ptr<MeshHandler> MeshHandlerPtr;

  class SurfaceMesh;
  typedef std::shared_ptr<SurfaceMesh> SurfaceMeshPtr;

  class LogicalVolume;
  typedef std::shared_ptr<LogicalVolume> LogicalVolumePtr;

  class FieldFunctionInterpolation;
  typedef FieldFunctionInterpolation FFInterp;
  typedef std::shared_ptr<FFInterp> FFInterpPtr;

  class UnpartitionedMesh;
  typedef std::shared_ptr<UnpartitionedMesh> UnpartitionedMeshPtr;
  typedef UnpartitionedMeshPtr UnpartMeshPtr;
}//namespace chi_mesh

namespace chi_physics
{
  class Solver;
  class Material;
  class TransportCrossSections;
  class FieldFunction;

  typedef std::shared_ptr<Solver>                 SolverPtr;
  typedef std::shared_ptr<Material>               MaterialPtr;
  typedef std::shared_ptr<TransportCrossSections> TransportCrossSectionsPtr;
  typedef std::shared_ptr<FieldFunction>          FieldFunctionPtr;
}//namespace chi_physics

namespace chi_math
{
  class Quadrature;
  class AngularQuadrature;

  typedef std::shared_ptr<Quadrature> QuadraturePtr;
  typedef std::shared_ptr<AngularQuadrature> AngularQuadraturePtr;

  class UnknownManager;
}//namespace chi_math

namespace chi_objects
{
  class MPI_Info;
  class ChiTimer;
  class ChiConsole;
  class ChiLog;
}//namespace chi_objects

//###################################################################
/**General utilities in ChiTech*/
class chi
{
public:
  static chi_objects::MPI_Info&   mpi;
  static chi_objects::ChiTimer    program_timer;
  static chi_objects::ChiConsole& console;
  static chi_objects::ChiLog&     log;

  static std::vector<chi_mesh::MeshHandlerPtr>   meshhandler_stack;
  static int         current_mesh_handler;

  static std::vector<chi_mesh::SurfaceMeshPtr>   surface_mesh_stack;
  static std::vector<chi_mesh::LogicalVolumePtr> logicvolume_stack;
  static std::vector<chi_mesh::FFInterpPtr>      field_func_interpolation_stack;
  static std::vector<chi_mesh::UnpartMeshPtr>    unpartitionedmesh_stack;

  static std::vector<chi_physics::SolverPtr>                 solver_stack;
  static std::vector<chi_physics::MaterialPtr>               material_stack;
  static std::vector<chi_physics::TransportCrossSectionsPtr> trnsprt_xs_stack;
  static std::vector<chi_physics::FieldFunctionPtr>          fieldfunc_stack;

  static std::vector<chi_math::QuadraturePtr>        quadrature_stack;
  static std::vector<chi_math::AngularQuadraturePtr> angular_quadrature_stack;



  //#######################################################
  /**Data block for run-time quantities.*/
  class run_time
  {
  public:
    static bool        termination_posted;
    static std::string input_file_name;
    static bool        sim_option_interactive;
    static bool        allow_petsc_error_handler;

  private:
    friend class chi;
    static void ParseArguments(int argc, char** argv);
    static int  InitPetSc(int argc, char** argv);
  public:

  public:
    run_time() = delete;                          //Deleted constructor
    run_time(const run_time&) = delete;           //Deleted copy constructor
    run_time operator=(const run_time&) = delete; //Deleted assigment operator
  };

  //#######################################################
  /**Customized exceptions.*/
  class RecoverableException : public std::runtime_error
  {
  public:
    explicit RecoverableException(const char* message) :
      std::runtime_error(std::string("RecoverableException: ") +
                         std::string(message)){}
    explicit RecoverableException(const std::string& message) :
      std::runtime_error(std::string("RecoverableException: ") +
                         message){}
     RecoverableException(const std::string& prefix,
                          const std::string& message) :
      std::runtime_error(prefix + message){}

    ~RecoverableException() noexcept override = default;
  };

public:
  chi() = delete;                     //Deleted constructor
  chi(const chi&) = delete;           //Deleted copy constructor
  chi operator=(const chi&) = delete; //Deleted assigment operator

public:
  static int  RunInteractive(int argc, char** argv);
  static int  RunBatch(int argc, char** argv);
  static int  Initialize(int argc, char** argv);
  static void Finalize();
  static void Exit(int error_code);

public:
  /**Attempts to retrieve an object of base-type `shared_ptr<T>` at the given
   * handle. It then attempts to cast it to type `shared_ptr<R>` and, if
   * successful, will return a reference of type R&.
   * \n
   * \n
   * Example usage:
   *
   * \code
   * const auto& surf_mesh = chi::GetStackItem<chi_mesh::SurfaceMesh>(
        chi::surface_mesh_stack, surface_hndl);
   * \endcode
   * */
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
      throw std::out_of_range("chi::GetStackItem: Invalid handle used. "
                             "Calling function: " + calling_function_name);
    }
  }

  /**Attempts to object of type `shared_ptr<T>` at the given
   * handle.
   * \n
   * \n
   * Example usage:
   *
   * \code
   * auto surf_mesh_ptr = chi::GetStackItemPtr<chi_mesh::SurfaceMesh>(
      chi::surface_mesh_stack, surf_mesh_hndle, fname);
   * \endcode
   * */
  template<class T>
  static std::shared_ptr<T>&
    GetStackItemPtr(std::vector<std::shared_ptr<T>>& stack,
                    const size_t handle,
                    const std::string& calling_function_name="Unknown")
  {
    try
    {
      std::shared_ptr<T>& item = stack.at(handle);
      return item;
    }
    catch (const std::out_of_range& oor)
    {
      throw std::out_of_range("chi::GetStackItem: Invalid handle used. "
                              "Calling function: " + calling_function_name);
    }
  }
};

#endif