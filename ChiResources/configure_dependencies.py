import sys
import os
import subprocess
import errno

print("########## Chi-Tech Dependency installer ##########")

####################################### Setting install dir
cwd = os.getcwd()
install_dir = f"{cwd}/chi-dependencies"
log_str = f">>{install_dir}/log.txt"
if len(sys.argv) == 1:
    print("Install directory not specified,", end="")
    print(f" defaulting to current directory:\n  \"{install_dir}\"")
else:
    cwd = sys.argv[1]
    if os.path.isdir(cwd):
        os.chdir(cwd)
        cwd = os.getcwd()
        install_dir = cwd
        print(f"Install directory set to \"{install_dir}\"")
    else:
        print("!!******!! Error !!*****!!: ", end="")
        print("Requested install path \"" + cwd + "\"", end="")
        print(" does not exist")
        sys.exit(1)
if not os.path.exists(install_dir):
    os.mkdir(install_dir)

####################################### Parsing versions
versions = {"readline": "8.0", "ncurses": "6.1", "lua": "5.3.5",
            "petsc": "3.12.5", "vtk": "8.2.0"}
packages = list(versions.keys())

for arg in sys.argv:
    for package in packages:
        req_version = versions[package]
        if package in arg:
            if "=" not in arg:
                print("To specify a version use \"<package>=<version>\".")
                print(f" Defaulting {package.upper()} "
                      f"to version {req_version}.")
                continue

            version = arg.split("=")[1]
            if version.count(".") != req_version.count("."):
                msg = f"Versions for {package.upper()} are specified via X"
                for _ in range(version.count(".")):
                    msg += ".X"
                print(msg)
                print(f" Defaulting {package.upper()} "
                      f"to version {req_version}.")
                continue

            req_version_vals = req_version.split(".")
            version_vals = version.split(".")
            zipped = zip(req_version_vals, version_vals)
            for i, (req_val, val) in enumerate(zipped):
                if int(val) < int(req_val):
                    print(f"Minimum requuired version for "
                          f"{package.upper()} not satisfied.")
                    print(f" Defaulting {package.upper()} "
                          f"to version {req_version}.")
                    continue
                if int(val) > int(req_val):
                    break

            versions[package] = version

readline_install = f"{install_dir}/READLINE/readline-{versions['ncurses']}/build"
ncurses_install = f"{install_dir}/NCURSES/ncurses-{versions['ncurses']}/build"
lua_install = f"{install_dir}/LUA/lua-{versions['lua']}/install"
petsc_install = f"{install_dir}/PETSC/petsc-{versions['petsc']}-install"
vtk_install = f"{install_dir}/VTK/VTK-{versions['vtk']}-install"

log_file = open(f"{install_dir}/log.txt", "w+")
roots_file = open(f"{install_dir}/configure_deproots.sh", "w+")


#######################################
# Prints a value to cout using system
def sprint(text):
    sys.stdout.write(text)


#######################################
# Prints a value to cout using system
# with prepended escape sequence
def sprint_esc(text):
    sys.stdout.write('\x1b' + text)


#######################################
# Runs a subprocess
def ExecSub(command, log, env_vars=None):
    success = True
    output = ""
    error = b"No Error"

    result = subprocess.Popen(command,
                              stdout=log,
                              stderr=subprocess.PIPE,
                              shell=True,
                              env=env_vars)
    output, error = result.communicate()

    if result.returncode != 0:
        success = False

    return success, error.decode("utf8")


#################################### Check for C Compiler
# Check for gcc and g++
def CheckForCCompilers():
    result = True
    print("Looking for c compiler...", end="")
    success, err = ExecSub("gcc --version", log_file)

    if success:
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
    print("Looking for c++ compiler...", end="")
    success, err = ExecSub("g++ --version", log_file)

    if success:
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
    print("Looking for fortran compiler...", end="")
    success, err = ExecSub("gfortran --version", log_file)

    if success:
        print("Success")
    else:
        print("Failed")
        print(err)
        result = False

    return result


####################################### Check dependency directory
def CheckDependencyDir():
    if not os.path.isdir(install_dir):
        print("Dependency directory not found")
    else:
        os.chdir(install_dir)


####################################### Get package
def DownloadPackage(url, pkg, ver, upper=False):
    if not os.path.exists(f"{install_dir}/{pkg.upper()}"):
        os.mkdir(f"{pkg.upper()}")
    os.chdir(f"{pkg.upper()}")

    pkg_ = pkg.upper() if upper else pkg
    if not os.path.exists(f"{os.getcwd()}/{pkg_}-{ver}.tar.gz"):
        print(f"Downloading {pkg_.upper()} {ver} to \"{os.getcwd()}\"")
        success, err = ExecSub(f"wget {url}", log_file)

        if upper:
            item = f"{pkg_}-{ver}.tar.gz"
            success, err = ExecSub(f"mv {item} {item.lower()}", log_file)


####################################### Get package
def ExtractPackage(pkg, ver):
    success, err = ExecSub(f"tar -zxf {pkg}-{ver}.tar.gz", log_file)

    os.chdir(f"{pkg}-{ver}")


####################################### Install readline
# readline_url = f"ftp://ftp.gnu.org/gnu/readline/readline-{versions['readline']}.tar.gz"
# readline_install = f"{install_dir}/READLINE/readline-{versions['readline']}/build"


def InstallReadline():
    pkg, ver = 'readline', versions['readline']
    url = f"ftp://ftp.gnu.org/gnu/{pkg}/{pkg}-{ver}.tar.gz"
    DownloadPackage(url, pkg, ver)

    # Check if it is installed already
    if not os.path.exists(f"{readline_install}/lib/libreadline.a"):
        ExtractPackage(pkg, ver)

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")

        env_vars = os.environ.copy()
        command = f"./configure --prefix={readline_install}"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)

        command = "make -j8"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)

        command = "make install"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)
    else:
        print(f"{pkg.capitalize()} already installed")

    os.chdir(install_dir)


####################################### Install ncurses
def InstallNcurses():
    pkg, ver = 'ncurses', versions['ncurses']
    url = f"https://invisible-mirror.net/archives/{pkg}/{pkg}-{ver}.tar.gz"
    DownloadPackage(url, pkg, ver)

    # Check if it is installed already
    if not os.path.exists(f"{ncurses_install}/lib/libncurses.a"):
        ExtractPackage(pkg, ver)

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")

        env_vars = os.environ.copy()
        command = f"./configure --prefix={ncurses_install}"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)

        command = "make -j8"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)

        command = "make install"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)
    else:
        print(f"{pkg.capitalize()} already installed")

    os.chdir(install_dir)


####################################### Install lua
# lua_url = f"https://www.lua.org/ftp/lua-{versions['lua']}.tar.gz"
# lua_install = f"{install_dir}/LUA/lua-{versions['lua']}/install"


def InstallLua():
    pkg, ver = 'lua', versions['lua']
    url = f"https://www.lua.org/ftp/{pkg}-{ver}.tar.gz"
    DownloadPackage(url, pkg, ver)

    # Check if it is installed already
    if not os.path.exists(f"{lua_install}/lib/liblua.a"):
        ExtractPackage(pkg, ver)

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")

        env_vars = os.environ.copy()

        lib_path = env_vars.get("LIBRARY_PATH", "")
        lib_path = f"{lib_path}:{readline_install}/lib:" \
                   f"{ncurses_install}/lib"
        env_vars["LIBRARY_PATH"] = lib_path

        c_path = env_vars.get('CPATH', "")
        c_path = f"{c_path}:{readline_install}/include"
        env_vars["CPATH"] = c_path

        os_tag = "linux"
        if "Darwin" in os.uname():
            os_tag = "macosx"

        command = f"make {os_tag} MYLIBS=-lncurses -j8"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)

        command = "make local"
        success, err = ExecSub(command, log_file, env_vars)
        if not success:
            print(command, err)
    else:
        print(f"{pkg.upper()} already installed")

    os.chdir(install_dir)


####################################### Install PETSc
# petsc_url = f"https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-{versions['petsc']}.tar.gz"
# petsc_install = f"{install_dir}/PETSC/petsc-{versions['petsc']}-install"


def InstallPETSc():
    pkg, ver = 'petsc', versions['petsc']
    url = f"https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/{pkg}-{ver}.tar.gz"
    DownloadPackage(url, pkg, ver)

    # Check if it is installed already
    if not os.path.exists(f"{petsc_install}/include/petsc"):
        ExtractPackage(pkg, ver)

        print(f"Configuring {pkg} {ver} to \"{os.getcwd()}\"")

        env_vars = os.environ.copy()
        env_vars["PETSC_DIR"] = f"{install_dir}/{pkg.upper()}/{pkg}-{ver}"

        exstring = f"./configure " \
                   f"--prefix={petsc_install} " \
                   f"--download-hypre=1 " \
                   f"--with-ssl=0 " \
                   f"--with-debugging=0 " \
                   f"--with-pic=1 " \
                   f"--with-shared-libraries=1 " \
                   f"--download-fblaslapack=1 " \
                   f"--download-metis=1 " \
                   f"--download-parmetis=1 " \
                   f"--download-superlu_dist=1 " \
                   f"--with-cxx-dialect=C++11 " \
                   f"--with-64-bit-indices " \
                   f"CFLAGS='-fPIC -fopenmp' " \
                   f"CXXFLAGS='-fPIC -fopenmp' " \
                   f"FFLAGS='-fPIC -fopenmp' " \
                   f"FCFLAGS='-fPIC -fopenmp' " \
                   f"F90FLAGS='-fPIC -fopenmp' " \
                   f"F77FLAGS='-fPIC -fopenmp' " \
                   f"COPTFLAGS='-O3 -march=native -mtune=native' " \
                   f"CXXOPTFLAGS='-O3 -march=native -mtune=native' " \
                   f"FOPTFLAGS='-O3 -march=native -mtune=native' " \
                   f"PETSC_DIR={env_vars['PETSC_DIR']}"

        print(exstring)
        success, err = ExecSub(exstring, log_file, env_vars)
        success, err = ExecSub("make all", log_file, env_vars)
        success, err = ExecSub("make install", log_file, env_vars)
    else:
        print(f"{pkg.upper()} already installed")

    os.chdir(install_dir)


####################################### Install VTK
# vtk_dir = ".".join(versions['vtk'].split(".")[:2])
# vtk_url = f"https://www.vtk.org/files/release/{vtk_dir}/VTK-{versions['vtk']}.tar.gz"
# vtk_install = f"{install_dir}/VTK/VTK-{versions['vtk']}-install"


def InstallVTK():
    pkg, ver = 'vtk', versions['vtk']
    vtk_dir = ".".join(ver.split(".")[:2])
    url = f"https://www.vtk.org/files/release/{vtk_dir}/{pkg.upper()}-{ver}.tar.gz"
    DownloadPackage(url, pkg, ver, upper=True)

    # Check if vtk is installed already
    if not os.path.exists(f"{vtk_install}/include"):
        ExtractPackage(pkg, ver)

        print(f"Configuring {pkg.upper()} {ver} to \"{os.getcwd()}\"")
        os.mkdir("build")
        os.chdir("build")

        env_vars = os.environ.copy()
        success, err = ExecSub(f"cmake -DCMAKE_INSTALL_PREFIX={vtk_install} "
                               f"-DBUILD_SHARED_LIBS:BOOL=ON "
                               f"-DVTK_Group_MPI:BOOL=ON "
                               f"-DVTK_GROUP_ENABLE_Qt=NO "
                               f"-DVTK_GROUP_ENABLE_Rendering=NO "
                               f"-DVTK_GROUP_ENABLE_Imaging=NO "
                               f"-DVTK_GROUP_ENABLE_StandAlone=WANT "
                               f"-DVTK_GROUP_ENABLE_Web=NO "
                               f"-DVTK_BUILD_TESTING:BOOL=OFF "
                               f"-DCMAKE_BUILD_TYPE=Release "
                               f"-DCMAKE_CXX_FLAGS=-std=c++11 "
                               f" ../", log_file, env_vars)
        success, err = ExecSub("make -j8", log_file, env_vars)
        success, err = ExecSub("make install", log_file, env_vars)
    else:
        print("VTK already installed")


success_c = CheckForCCompilers()
if not success_c:
    print("***** c compiler not working")
    exit(errno.EPERM)
success_cpp = CheckForCPPCompilers()
if not success_cpp:
    print("***** c++ compiler not working")
    exit(errno.EPERM)
success_f = CheckFortranCompilers()
if not success_f:
    print("***** Fortran compiler not working")
    exit(errno.EPERM)

CheckDependencyDir()
InstallReadline()
InstallNcurses()
InstallLua()

# InstallBoost()
InstallPETSc()
InstallVTK()

roots_file.write('export BASE_PATH="' + install_dir + '"\n')
roots_file.write('\n')
roots_file.write(f'export LUA_ROOT="{lua_install}"\n')
roots_file.write(f'export PETSC_ROOT="{petsc_install}"\n')
roots_file.write(f'export VTK_DIR="{vtk_install}"\n')
roots_file.write(f'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"{vtk_install}/lib"\n')
roots_file.write('echo "Enviroment set for compiling. If recompiling changed sources execute"\n')
roots_file.write('echo "     ./configure clean"\n')
roots_file.write('echo " "\n')
roots_file.write('echo "Otherwise just execute:"\n')
roots_file.write('echo "     ./configure"\n')

print(os.getcwd(), os.listdir())
ExecSub("chmod u+x configure_deproots.sh", log_file)

log_file.close()
roots_file.close()

print("########## Chi-Tech Dependency install complete ##########")
print(f"Now execute: \n     $. {install_dir}/configure_deproots.sh\n")
