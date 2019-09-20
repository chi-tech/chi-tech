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

#### Step 2 - Install Boost

Download and extract **Boost** from https://www.boost.org.  **Boost** v1.63.0+ is
supported.  Configure and install as follows:
```console
    $ ./bootstrap.sh --prefix=/Path/to/Boost
    $ ./b2 install
```

Once the installation is complete, set the BOOST_ROOT environment variable to the
**BOOST** install location:
```console
    $ export BOOST_ROOT=/Path/to/Boost
```
#### Step 4 - Install PETSc

Download and extract **PETSc** from https://www.mcs.anl.gov/petsc.
**PETSc** v3.9.4 is recommended for performance reasons.  Configure and install
as follows:
```console
    $ ./configure --prefix=/Path/to/PETSc \
    --with-fortran-bindings=0 \
    --download-hypre=1 \
    --with-ssl=0 \
    --with-debugging=0 \
    --with-pic=1 \
    --with-shared-libraries=1 \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack=1 \
    --with-cxx-dialect=C++11 \
    CFLAGS='-fPIC' \
    CXXFLAGS='-fPIC' \
    FCFLAGS='-fPIC' \
    COPTFLAGS='-O3' \
    CXXOPTFLAGS='-O3' \
    FOPTFLAGS='-O3'
```
If the configuration fails consult the **PETSc** user documentation.  Upon
completion of the configure step, **PETSc** will provide the make command
that you should use. Execute this command.  After a successful make, **PETSc**
will provide the command used to install the header files and libraries.
Execute this command.

To test whether **PETSc** has been installed correctly execute:
```console
    $ make test
```
Finally, set the PETSC_ROOT environment variable to the **PETSc** install location:
```console
    $ export PETSC_ROOT=/Path/to/PETSC
```
#### Step 5 - Install the Visualization Tool Kit

Download and extract **VTK** from https://vtk.org.  **VTK** v8.2.0+ is recommended.
Configure and install as follows:
```console
    $ mkdir vtk-build
    $ cd vtk-build
    $ cmake -DCMAKE_INSTALL_PREFIX=/Path/to/VTK \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DVTK_Group_MPI:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
    /Path/to/VTK/source

    $ make
    $ make install
```

Set the VTK_DIR environment variable to the **VTK** install location:
```console
    $ export VTK_DIR=/Path/to/VTK
```

#### Step 6 - Install Eigen

Download and extract **Eigen** from https://eigen.tuxfamily.org.  **Eigen** v3.3.7+
is recommended.  **Eigen** is a header only library, and no further installation is
required.

Set the EIGEN_ROOT environment variable to the **Eigen** install location:
```console
    $ export EIGEN_ROOT=/Path/to/Eigen
```

#### Step 7 - Install Random123

Download and extract **Random123** from https://www.deshawresearch.com/resources_random123.html.
**Radom123** v1.13.2+ is recommended.  **Random123** is a header only library,
and no further installation is required.

Set the RANDOM123_ROOT environment variable to the **Random123** install location:
```console
    $ export RANDOM123_ROOT=/Path/to/Random123
```

#### Step 8 - Install Triangle

Download and extract **Triangle** from https://www.cs.cmu.edu/~quake/triangle.html.
**Triangle** v1.6 is recommended.  Before installing **Triangle** edit the makefile
and modify the line that reads:
```console
    CSWITCHES = -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib
```
to:
```console
    CSWITCHES = -O -I/usr/X11R6/include -L/usr/X11R6/lib
```

Install **Triangle** as follows:
```console
    $ make
    $ make trilibrary
```

Set the TRIANGLE_ROOT environment variable to the **Triangle** install location:
```console
    $ export TRIANGLE_ROOT=/Path/to/Triangle
```

#### Step 9 - Install Lua

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

#### Step 10 - Configure and Build ChiTech

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

#### Step 11 - Generate the documentation

You can either access the documentation online [here](https://chi-tech.github.io),
or generate it locally.

To generate the documentation locally, Google and follow the online instructions
for installing installing doxygen and lua then run:
```console
    $ cd CHI_DOC/
    $ ./YReGenerateDocumentation.sh
```
Open the documentation index file at CHI_DOC/HTMLdocs/html/index.html.
