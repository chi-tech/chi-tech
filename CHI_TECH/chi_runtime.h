#ifndef chi_runtime_h
#define chi_runtime_h

#include <vector>
#include <string>

namespace chi_mesh
{
  class MeshHandler;
}

/**General utilities in ChiTech*/
class ChiTech
{
public:
  static std::vector<chi_mesh::MeshHandler*>  meshhandler_stack;
  static int         current_mesh_handler;
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

  ChiTech() = delete;
};

#endif