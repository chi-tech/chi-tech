import sys
import os
import subprocess
import errno

print("########## Chi-Tech Dependency installer ##########")

####################################### Setting install dir
cwd = os.getcwd()
install_dir = cwd + "/chi-dependencies"
log_str = ">>" + install_dir + "log.txt"
if (len(sys.argv) == 1):
  print("Install directory not specified,",end="")
  print(" defaulting to current directory:\n  \"" + install_dir + "\"")
else:
  cwd = sys.argv[1]
  if (os.path.isdir(cwd)):
    os.chdir(cwd)
    cwd = os.getcwd()
    install_dir = cwd + "/"
    print("Install directory set to \"" + install_dir + "\"")
  else:
    print("!!******!! Error !!*****!!: ", end="")
    print("Requested install path \"" + cwd + "\"", end="")
    print(" does not exist")
    sys.exit(1)
if (not os.path.exists(install_dir)):
  os.mkdir(install_dir)

log_file = open(install_dir + "/log.txt","w+")
roots_file = open(install_dir + "/configure_deproots.sh","w+")

#######################################
# Prints a value to cout using system
def sprint(text):
  sys.stdout.write(text)

#######################################
# Prints a value to cout using system
# with prepended escape sequence
def sprint_esc(text):
  sys.stdout.write('\x1b'+text)

#######################################
# Runs a subprocess
def ExecSub(command,log,env_vars=None):
  success = True
  output = ""
  error = b"No Error"

  result = subprocess.Popen(command,
                            stdout=log,
                            stderr=subprocess.PIPE,
                            shell=True,
                            env=env_vars)
  output, error = result.communicate()

  if (result.returncode != 0):
    success = False

  return success,error.decode("utf8")

#################################### Check for C Compiler
# Check for gcc and g++
def CheckForCCompilers():
  result = True
  print("Looking for c compiler...",end="")
  success,err = ExecSub("gcc --version",log_file)

  if (success):
    print("Success")
  else:
    print("Failed")
    print(err)
    result = False

  return result

#################################### Check for C++ compiler
# Check for gcc and g++
def CheckForCPPCompilers():
  result = True
  print("Looking for c++ compiler...",end="")
  success,err = ExecSub("g++ --version",log_file)

  if (success):
    print("Success")
  else:
    print("Failed")
    print(err)
    result = False

  return result

#################################### Check for fortran compiler
# Check for gcc and g++
def CheckFortranCompilers():
  result = True
  print("Looking for fortran compiler...",end="")
  success,err = ExecSub("gfortran --version",log_file)

  if (success):
    print("Success")
  else:
    print("Failed")
    print(err)
    result = False

  return result

####################################### Check dependency directory
def CheckDependencyDir():
  if (not os.path.isdir(install_dir)):
    print("Dependency directory not found")
  else:
    os.chdir(install_dir)

####################################### Install Eigen
eigen_url = "http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz"
def InstallEigen():
  if (not os.path.exists(install_dir + "/EIGEN")):
    os.mkdir("EIGEN")
  os.chdir("EIGEN")

  if (not os.path.exists(install_dir + "/EIGEN/3.3.7.tar.gz")):
    print("Downloading Eigen 3.3.7 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + eigen_url,log_file)

  #Check if eigen is installed already
  installed = os.path.exists(install_dir + "/EIGEN/eigen-3.3.7/Eigen")
  if (not installed):
    print("Configuring Eigen 3.3.7 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf 3.3.7.tar.gz",log_file)
    success,err = ExecSub("mv eigen-* eigen-3.3.7",log_file)
  else:
    print("EIGEN already installed")

  os.chdir(install_dir)

####################################### Install Random-123
r123_url = "https://www.deshawresearch.com/downloads/download_random123.cgi/Random123-1.13.2.tar.gz"
def InstallRandom123():
  if (not os.path.exists(install_dir + "/RANDOM123")):
    os.mkdir("RANDOM123")
  os.chdir("RANDOM123")

  if (not os.path.exists(install_dir + "/RANDOM123/Random123-1.13.2.tar.gz")):
    print("Downloading Random123-1.13.2 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + r123_url,log_file)

  #Check if it is installed already
  installed = os.path.exists(install_dir + "/RANDOM123/Random123-1.13.2/include")
  if (not installed):
    print("Configuring Random123-1.13.2 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf Random123-1.13.2.tar.gz",log_file)
  else:
    print("RANDOM123 already installed")

  os.chdir(install_dir)

####################################### Install triangle
tri_url = "http://www.netlib.org/voronoi/triangle.zip"
def InstallTriangle():
  if (not os.path.exists(install_dir + "/TRIANGLE")):
    os.mkdir("TRIANGLE")
  os.chdir("TRIANGLE")

  if (not os.path.exists(install_dir + "/TRIANGLE/triangle.zip")):
    print("Downloading Triangle to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + tri_url,log_file)

  #Check if it is installed already
  installed = os.path.exists(install_dir + "/TRIANGLE/triangle.o")
  if (not installed):
    print("Configuring Triangle to \"" + os.getcwd() + "\"")
    success,err = ExecSub("unzip triangle.zip",log_file)

    os_tag = "LINUX"
    if ("Darwin" in os.uname() ):
      os_tag = "APPLE"

    env_vars=os.environ.copy()
    command = "gcc -O -D" + os_tag + " -I/usr/X11R6/include" + \
              " -L/usr/X11R6/lib -o ./triangle " + \
              " ./triangle.c -lm"
    success,err = ExecSub(command,log_file,env_vars)
    print(command,err)

    command = "gcc -O -D" + os_tag + " -I/usr/X11R6/include " + \
              "-L/usr/X11R6/lib -DTRILIBRARY -c -o " + \
              "./triangle.o ./triangle.c"
    success,err = ExecSub(command,log_file,env_vars)
    print(command,err)
  else:
    print("Triangle already installed")

  os.chdir(install_dir)

####################################### Install readline
readline_url = "ftp://ftp.gnu.org/gnu/readline/readline-8.0.tar.gz"
def InstallReadline():
  if (not os.path.exists(install_dir + "/READLINE")):
    os.mkdir("READLINE")
  os.chdir("READLINE")

  if (not os.path.exists(install_dir + "/READLINE/readline-8.0.tar.gz")):
    print("Downloading Readline 8.0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + readline_url,log_file)

  #Check if it is installed already
  installed = os.path.exists(install_dir + "/READLINE/readline-8.0"
                                           "/build/lib/libreadline.a")
  if (not installed):
    print("Configuring Readline 8.0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf readline-8.0.tar.gz",log_file)

    os.chdir("readline-8.0")
    env_vars=os.environ.copy()
    command = "./configure --prefix=" + \
              install_dir + \
              "/READLINE/readline-8.0/build"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)

    command = "make -j8"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)

    command = "make install"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)
  else:
    print("Readline already installed")

  os.chdir(install_dir)

####################################### Install ncurses
ncurses_url = "https://invisible-mirror.net/archives/ncurses/ncurses-6.1.tar.gz"
lua_url = "https://www.lua.org/ftp/lua-5.3.5.tar.gz"
def Install_ncurses():
  if (not os.path.exists(install_dir + "/NCURSES")):
    os.mkdir("NCURSES")
  os.chdir("NCURSES")

  if (not os.path.exists(install_dir + "/NCURSES/ncurses-6.1.tar.gz")):
    print("Downloading NCurses 6.1 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + ncurses_url,log_file)

  #Check if it is installed already
  installed = os.path.exists(install_dir + "/NCURSES/ncurses-6.1"
                                           "/build/lib/libncurses.a")
  if (not installed):
    print("Configuring ncurses 6.1 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf ncurses-6.1.tar.gz",log_file)

    os.chdir("ncurses-6.1")
    env_vars=os.environ.copy()
    command = "./configure --prefix=" + \
              install_dir + \
              "/NCURSES/ncurses-6.1/build"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)

    command = "make -j8"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)

    command = "make install"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)
  else:
    print("Ncurses already installed")

  os.chdir(install_dir)

####################################### Install lua
lua_url = "https://www.lua.org/ftp/lua-5.3.5.tar.gz"
def InstallLua():
  if (not os.path.exists(install_dir + "/LUA")):
    os.mkdir("LUA")
  os.chdir("LUA")

  if (not os.path.exists(install_dir + "/LUA/lua-5.3.5.tar.gz")):
    print("Downloading Lua 5.3.5 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + lua_url,log_file)

  #Check if it is installed already
  installed = os.path.exists(install_dir + "/LUA/lua-5.3.5"
                                           "/install/lib/liblua.a")
  if (not installed):
    print("Configuring LUA 5.3.5 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf lua-5.3.5.tar.gz",log_file)

    os.chdir("lua-5.3.5")
    env_vars=os.environ.copy()

    lib_path = env_vars.get("LIBRARY_PATH","")
    lib_path = lib_path + ":" + install_dir + "/READLINE/readline-8.0/build/lib"
    lib_path = lib_path + ":" + install_dir + "/NCURSES/ncurses-6.1/build/lib"
    env_vars["LIBRARY_PATH"] = lib_path

    c_path =  env_vars.get('CPATH',"")
    c_path = c_path + ":" + install_dir + "/READLINE/readline-8.0/build/include"
    env_vars["CPATH"] = c_path

    os_tag = "linux"
    if ("Darwin" in os.uname() ):
      os_tag = "macosx"

    command = "make " + os_tag + " MYLIBS=-lncurses -j8"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)
    print(command,err)

    command = "make local"
    success,err = ExecSub(command,log_file,env_vars)
    if (not success):
      print(command,err)
  else:
    print("Lua already installed")

  os.chdir(install_dir)


####################################### Install boost
boost_url = "https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.gz"
def InstallBoost():
  if (not os.path.exists(install_dir + "/BOOST")):
    os.mkdir("BOOST")
  os.chdir("BOOST")

  if (not os.path.exists(install_dir + "/BOOST/boost_1_71_0.tar.gz")):
    print("Downloading Boost 1.71_0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + boost_url,log_file)
    sprint_esc("[1A")
    sprint_esc("[2K")

  #Check if boost is installed already
  installed = os.path.exists(install_dir + "/BOOST/boost_1_71_0/include")
  if (not installed):
    print("Configuring Boost 1.71_0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf boost_1_71_0.tar.gz",log_file)
    os.chdir("boost_1_71_0/")
    success,err = ExecSub("./bootstrap.sh --with-toolset=gcc",log_file)
    success,err = ExecSub("./b2 headers",log_file)
    os.mkdir("include")
    success,err = ExecSub("cp -r boost include/",log_file)
  else:
    print("BOOST already installed")

  os.chdir(install_dir)

####################################### Install PETSc
petsc_url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz"
def InstallPETSc():
  if (not os.path.exists(install_dir + "/PETSc")):
    os.mkdir("PETSc")
  os.chdir("PETSc")

  if (not os.path.exists(install_dir + "/PETSc/petsc-3.9.4.tar.gz")):
    print("Downloading PETSc 3.9.4 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + petsc_url,log_file)


  env_vars=os.environ.copy()
  env_vars["PETSC_DIR"] = install_dir + "/PETSc/petsc-3.9.4"

  #Check if petsc is installed already
  installed = os.path.exists(install_dir + \
                             "/PETSc/petsc-3.9.4/install/include/petsc")
  if (not installed):
    print("Configuring PETSc 3.9.4 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf petsc-3.9.4.tar.gz",log_file)
    os.chdir("petsc-3.9.4/")
    exstring = "./configure " \
        "--prefix=" + install_dir + "/PETSc/petsc-3.9.4/install " \
        "--download-hypre=git://https://github.com/hypre-space/hypre.git " \
        "--with-ssl=0 " \
        "--with-debugging=0 " \
        "--with-pic=1 " \
        "--with-shared-libraries=1 " \
        "--download-fblaslapack=1 " \
        "--download-metis=1 " \
        "--download-parmetis=1 " \
        "--download-superlu_dist=1 " \
        "--with-cxx-dialect=C++11 " \
        "CFLAGS='-fPIC -fopenmp' " \
        "CXXFLAGS='-fPIC -fopenmp' " \
        "FFLAGS='-fPIC -fopenmp' " \
        "FCFLAGS='-fPIC -fopenmp' " \
        "F90FLAGS='-fPIC -fopenmp' " \
        "F77FLAGS='-fPIC -fopenmp' " \
        "COPTFLAGS='-O3 -march=native -mtune=native' " \
        "CXXOPTFLAGS='-O3 -march=native -mtune=native' " \
        "FOPTFLAGS='-O3 -march=native -mtune=native' " \
        "PETSC_DIR=" + install_dir + "/PETSc/petsc-3.9.4/"

    print(exstring)
    success,err = ExecSub(exstring,log_file,env_vars)
    success,err = ExecSub("make all",log_file,env_vars)
    success,err = ExecSub("make install",log_file,env_vars)
  else:
    print("PETSc already installed")

  os.chdir(install_dir)


####################################### Install VTK
vtk_url   = "https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz"
def InstallVTK():
  if (not os.path.exists(install_dir + "/VTK")):
    os.mkdir("VTK")
  os.chdir("VTK")

  if (not os.path.exists(install_dir + "/VTK/VTK-8.2.0.tar.gz")):
    print("Downloading VTK 8.2.0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("wget " + vtk_url,log_file)


  # Check if vtk is installed already
  installed = os.path.exists(install_dir + "/VTK/VTK-8.2.0/install/include")
  if (not installed):
    print("Configuring VTK 8.2.0 to \"" + os.getcwd() + "\"")
    success,err = ExecSub("tar -zxf VTK-8.2.0.tar.gz",log_file)
    os.chdir("VTK-8.2.0/")
    os.mkdir("build")
    os.chdir("build")

    env_vars=os.environ.copy()
    success,err = ExecSub("cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install " + \
                          "-DBUILD_SHARED_LIBS:BOOL=ON " + \
                          "-DVTK_Group_MPI:BOOL=ON " + \
                          "-DVTK_GROUP_ENABLE_Qt=NO " + \
                          "-DVTK_GROUP_ENABLE_Rendering=NO " + \
                          "-DVTK_GROUP_ENABLE_Imaging=NO " + \
                          "-DVTK_GROUP_ENABLE_StandAlone=WANT " + \
                          "-DVTK_GROUP_ENABLE_Web=NO " + \
                          "-DVTK_BUILD_TESTING:BOOL=OFF " + \
                          "-DCMAKE_BUILD_TYPE=Release " + \
                          "-DCMAKE_CXX_FLAGS=-std=c++11 " + \
                          " ../",log_file,env_vars)
    success,err = ExecSub("make -j8",log_file,env_vars)
    success,err = ExecSub("make install",log_file,env_vars)
  else:
    print("VTK already installed")


success = CheckForCCompilers()
if (not success):
  print("***** c compiler not working")
  exit(errno.EPERM)
success = CheckForCPPCompilers()
if (not success):
  print("***** c++ compiler not working")
  exit(errno.EPERM)
success = CheckFortranCompilers()
if (not success):
  print("***** Fortran compiler not working")
  exit(errno.EPERM)



CheckDependencyDir()
InstallEigen()
InstallRandom123()
InstallTriangle()
InstallReadline()
Install_ncurses()
InstallLua()

# InstallBoost()
InstallPETSc()
InstallVTK()



roots_file.write('export BASE_PATH="' + install_dir + '"\n')
roots_file.write('\n')
roots_file.write('export EIGEN_ROOT="$BASE_PATH/EIGEN/eigen-3.3.7"\n')
# roots_file.write('export RANDOM123_ROOT="$BASE_PATH/RANDOM123/Random123-1.13.2"\n')
# roots_file.write('export TRIANGLE_ROOT="/$BASE_PATH/TRIANGLE"\n')
roots_file.write('export LUA_ROOT="$BASE_PATH/LUA/lua-5.3.5/install"\n')
# roots_file.write('export BOOST_ROOT="$BASE_PATH/BOOST/boost_1_71_0"\n')
roots_file.write('export PETSC_ROOT="$BASE_PATH/PETSc/petsc-3.9.4/install"\n')
roots_file.write('export VTK_DIR="$BASE_PATH/VTK/VTK-8.2.0/install"\n')
roots_file.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$BASE_PATH/VTK/VTK-8.2.0/install/lib"\n')
roots_file.write('echo "Enviroment set for compiling. If recompiling changed sources execute"\n')
roots_file.write('echo "     ./configure clean"\n')
roots_file.write('echo " "\n')
roots_file.write('echo "Otherwise just execute:"\n')
roots_file.write('echo "     ./configure"\n')

ExecSub("chmod u+x configure_deproots.sh",log_file)

log_file.close()
roots_file.close()

print("########## Chi-Tech Dependency install complete ##########")
print("Now execute: \n     $. " + install_dir + "/configure_deproots.sh\n")
