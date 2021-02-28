## Compiling on macOS
___
### Compiling with the Apple Toolchain
___

The following instructions assume that all of the prerequisite software packages
will be installed by hand.  These packages may also be installed via **Brew** or
**MacPorts**.

#### Step 1 - Install MPI

Install either **OpenMPI** or **MPICH**.  **MPICH** is recommended for
better performance and stability.

**MPICH**: https://www.mpich.org/

**OpenMPI**: https://www.open-mpi.org/

Download and extract your preferred MPI package and install:
```console
    $ ./configure --prefix=/Path/to/MPI
    $ make
    $ make install
```

Once the installation is complete, add the MPI binaries directory to your path:
```console
    $ export PATH=/Path/to/MPI/bin:$PATH
```
<u>NOTE:</u> *You may want to permanently add this and the other `export ...`
commands below to your bash profile file, so that you don't have to execute these
commands every time you open a new terminal. This is typically done by adding
these commands to the file `.profile` or `.bash_profile` in your home directory.*

To check that the install was successful, execute one of the compiler wrappers:
```console
    $ mpicxx -v
    mpicxx for MPICH version 3.3
    Apple LLVM version 10.0.1 (clang-1001.0.46.4)
    Target: x86_64-apple-darwin18.7.0
    Thread model: posix
    InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin
```

#### Step 2 - Install PETSc

The current supported version is
[petsc version 3.12.5](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.12.5.tar.gz) .


Return to your *projects* folder (or whatever you chose to place stuff). Run
the following

```bash
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.12.5.tar.gz
tar -zxf petsc-3.12.5.tar.gz
cd petsc-3.12.5
```

Download and extract the archive file to a folder of your choice then navigate
into the directory containing the "configure" script and execute the following:

```bash
./configure  \
--prefix=$PWD/install  \
--download-hypre=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-shared-libraries=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--with-cxx-dialect=C++11  \
CFLAGS='-fPIC -fopenmp'  \
CXXFLAGS='-fPIC -fopenmp'  \
FFLAGS='-fPIC -fopenmp'  \
FCFLAGS='-fPIC -fopenmp'  \
F90FLAGS='-fPIC -fopenmp'  \
F77FLAGS='-fPIC -fopenmp'  \
COPTFLAGS='-O3 -march=native -mtune=native'  \
CXXOPTFLAGS='-O3 -march=native -mtune=native'  \
FOPTFLAGS='-O3 -march=native -mtune=native'  \
PETSC_DIR=$PWD
```

If the configuration fails then consult PETSc's user documentation.

<u>NOTE:</u> *The recommended version of PETSc requires Python 2.6+ for its configuration.
If you need to use Python 3.4+ (which is typically pre-installed on modern Linux systems as
`python3`), you will need to download PETSc 3.11 (current most recent version is 3.11.3)
and run the above configure script as `python3 ./configure ...`*

Upon completion of the configure step, PETSc will provide the make command
that you should use. Execute this command.

After a successful make PETSc should indicate the command used to install
the header files and libraries. Execute the associated command as well.

To test whether the system has been installed correctly execute:

```bash
make test
```

The final step is to generate the environment variable for the install
directory of PETSc:

```bash
export PETSC_ROOT=$PWD/install
```

Again, this is also something you'd like to add to your bash profile.

#### Step 3 - Install the Visualization Tool Kit

In your projects folder install VTK using the following commands:

```bash
mkdir VTK
cd VTK
wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar -zxf VTK-8.2.0.tar.gz
cd VTK-8.2.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install  + \
-DBUILD_SHARED_LIBS:BOOL=ON  + \
-DVTK_Group_MPI:BOOL=ON  + \
-DVTK_GROUP_ENABLE_Qt=NO  + \
-DVTK_GROUP_ENABLE_Rendering=NO  + \
-DVTK_GROUP_ENABLE_Imaging=NO  + \
-DVTK_GROUP_ENABLE_StandAlone=WANT  + \
-DVTK_GROUP_ENABLE_Web=NO  + \
-DVTK_BUILD_TESTING:BOOL=OFF  + \
-DCMAKE_BUILD_TYPE=Release  + \
-DCMAKE_CXX_FLAGS=-std=c++11  + \
 ../
```

This will take a while after which you need to execute:

```bash
make -j4 && make install
```

This will also take a while. Finally, export the *VTK_DIR* variable:

```bash
export VTK_DIR=$PWD/install
```

#### Step 4 - Install Lua

Download and extract **Lua** from https://www.lua.org.  v5.3.5+ is recommended.
Before installing **Lua** edit the Makefile and set INSTALL_TOP to your desired
install location.  Install **Lua** as follows:
```console
    $ make macosx
    $ make install
```
If the install complains about missing **readline** includes or libraries, it may
be necessary to install **readline** first.

Set the LUA_ROOT environment variable to the **Lua** install location:
```console
    $export LUA_ROOT=/Path/to/Lua
```

#### Step 5 - Configure and Build ChiTech

Clone the **ChiTech** repository:
```console
    $ git clone https://github.com/chi-tech/chi-tech
```

Run the configure script:
```console
    $ cd chi-tech
    $ ./configure.sh
```
The configure script will generate the CMake build scripts.

Once configure.sh has completed, build **ChiTech** with:
```console
    $ make
```

### Step 6 - Run regression tests

To check if the code compiled correctly execute the test scripts:

```bash
    $ python3 CHI_TEST/Z_Run_all.py
```


#### Step 7 - Generate the documentation

You can either access the documentation online [here](https://chi-tech.github.io),
or generate it locally.

To generate the documentation locally, Google and follow the online instructions
for installing installing doxygen and lua then run:
```console
    $ cd CHI_DOC/
    $ ./YReGenerateDocumentation.sh
```
Open the documentation index file at CHI_DOC/HTMLdocs/html/index.html.
