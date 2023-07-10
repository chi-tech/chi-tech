import argparse
import os
import pathlib
import shutil
import stat

print("**********")
posix_path = pathlib.Path(__file__).parent.resolve()
chi_tech_root = os.path.abspath(str(posix_path) + "/../")
print("ChiTech root: ", chi_tech_root)

cwd = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("module_name", type=str,
                    help="The name of the module")

args = parser.parse_args()

mod_full_name = args.module_name

[mod_path, mod_ext] = os.path.splitext(mod_full_name)
mod_name = os.path.basename(mod_path)

print("**********")
print("Creating module with name \"" + mod_name + "\"")
print("in folder \"" + mod_path + "\"")

# =============================================== Create folder structure
dirs = mod_path.split("/")
incremental_dir = "./"
for dur in dirs:
    incremental_dir += "/" + dur
    if not os.path.exists(incremental_dir):
        os.mkdir(incremental_dir)

if not os.path.exists(mod_path):
    os.mkdir(mod_path)

if not os.path.exists(mod_path + "/modules"):
    os.mkdir(mod_path + "/modules")

if not os.path.exists(mod_path + "/modules/" + mod_name):
    os.mkdir(mod_path + "/modules/" + mod_name)

if not os.path.exists(mod_path + "/tests"):
    os.mkdir(mod_path + "/tests")

# =============================================== Copy files
shutil.copyfile(chi_tech_root+"/configure.sh", mod_path+"/configure.sh")
st = os.stat(mod_path+"/configure.sh")
os.chmod(mod_path+"/configure.sh", st.st_mode | stat.S_IEXEC)

# =============================================== Create files
# test.lua
test_dot_lua = open(mod_path+"/tests/test.lua", "w")
test_dot_lua.write("chiLog(LOG_0, \"Hello World!\")\n")
test_dot_lua.close()

# main.cc
main_dot_cc = open(mod_path+"/main.cc", "w")
main_dot_cc.write('''\
#include "chi_runtime.h"

//######################################################### Program entry point
/** Program entry point.

\\param argc int    Number of arguments supplied.
\\param argv char** Array of strings representing each argument.

*/
int main(int argc, char* argv[])
{
  Chi::Initialize(argc,argv, MPI_COMM_WORLD);
  
  int error_code;
  if (Chi::run_time::sim_option_interactive_)
    error_code = Chi::RunInteractive(argc, argv);
  else
    error_code = Chi::RunBatch(argc, argv);

  Chi::Finalize();

  return error_code;
}
''')
main_dot_cc.close()

# CMakeLists.txt
cml_dot_txt = open(mod_path+"/CMakeLists.txt", "w")
cml_dot_txt.write('''\
cmake_minimum_required(VERSION 3.12)
\n''' +
f"set(TARGET {mod_name})\n" +
f"project({mod_name} LANGUAGES CXX)" + '''

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED CHI_TECH_DIR)
    if (NOT (DEFINED ENV{CHI_TECH_DIR}))
        message(FATAL_ERROR "***** CHI_TECH_DIR is not set *****")
    else()
        set(CHI_TECH_DIR "$ENV{CHI_TECH_DIR}")
    endif()
endif()
message(STATUS "CHI_TECH_DIR set to ${CHI_TECH_DIR}")

include("${CHI_TECH_DIR}/resources/CMakeMacros/Downstream.cmake")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

file (GLOB_RECURSE SOURCES "*.cc")
add_executable(${TARGET} "${SOURCES}")
if(UNIX AND NOT APPLE)
    target_link_libraries(${TARGET} ${CHI_LIBS} -Wl,--whole-archive ChiLib -Wl,--no-whole-archive )
elseif(APPLE)
    target_link_libraries(${TARGET} ${CHI_LIBS} -Wl,-all_load ChiLib )
endif()

file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build \n\n"
        "clean:\n\t$(MAKE) -C chi_build clean\n")
''')

cml_dot_txt.close()

