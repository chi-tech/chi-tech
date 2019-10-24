import sys
import time
import os
import subprocess
import multiprocessing

boost_url = "https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.gz"
petsc_url = "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz"

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
    install_dir = cwd + "/chi-dependencies"
    print("Install directory set to \"" + install_dir + "\"")
  else:
    print("!!******!! Error !!*****!!: ", end="")
    print("Requested install path \"" + cwd + "\"", end="")
    print(" does not exist")
    sys.exit(1)
if (not os.path.exists(install_dir)):
  os.mkdir("chi-dependencies")

log_file = open(install_dir + "/log.txt","w+")

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

####################################### Install boost
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
        "--download-hypre=1 " \
        "--with-ssl=0 " \
        "--with-debugging=0 " \
        "--with-pic=1 " \
        "--with-shared-libraries=1 " \
        "--download-fblaslapack=1 " \
        "--download-metis=1 " \
        "--download-parmetis=1 " \
        "--download-superlu_dist=1 " \
        "--with-cxx-dialect=C++11 " \
        "--download-mpich=1 "\
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
  else:
    print("PETSc already installed")

CheckForCCompilers()
CheckForCPPCompilers()
CheckFortranCompilers()
CheckDependencyDir()
InstallBoost()
InstallPETSc()

log_file.close()

########## FOR DEMO ################
if __name__ == "__main__":
  print("Mother..0%")
  for i in range(1,101):
      time.sleep(0.01)
      sprint_esc("[1A")
      sprint_esc("[2K")
      print("Mother.." + str(i) + "%")



####################################

