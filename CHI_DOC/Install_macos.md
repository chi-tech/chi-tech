## Compiling on macOS
___
### Compiling with the Apple Toolchain
___

The following instructions assume that all of the prerequisite software packages will be installed by hand.  These packages may also be installed via **Brew** or **MacPorts**.

#### Step 1 - Install MPI

Install either **OpenMPI** or **MPICH**.  **MPICH** is recommended for
better performance and stability.

**MPICH**: https://www.mpich.org/

**OpenMPI**: https://www.open-mpi.org/

Download and extract your preferred MPI package and install:

    $ ./configure --prefix=/Path/to/MPI
    $ make
    $ make install

Once the installation is complete, add the /Path/to/MPI/bin directory to your path:

    $ export PATH=/Path/to/MPI/bin:$PATH

<u>NOTE:</u> *You may want to permanently add this and the other `export ...` commands below to 
your bash profile file, so that you don't have to execute these commands every time you open a new 
terminal. This is typically done by adding these commands to the file `.profile` or `.bash_profile` 
in your home directory.*

To check that the install was successful, execute one of the compiler wrappers:

    $ mpicxx -v
    mpicxx for MPICH version 3.3
    Apple LLVM version 10.0.1 (clang-1001.0.46.4)
    Target: x86_64-apple-darwin18.7.0
    Thread model: posix
    InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin

#### Step 2 - Install Boost

Download and extract **Boost** from https://www.boost.org.  **Boost** v1.63.0+ is supported.  Configure and install as follows:

    $ ./bootstrap.sh --prefix=/Path/to/Boost
    $ ./b2 install

Once the installation is complete, set the BOOST_ROOT environment variable:

    $ export BOOST_ROOT=/Path/to/Boost

#### Step 4 - Install PETSc

Download and extract **PETSc** from https://www.mcs.anl.gov/petsc.  **PETSc** v3.9.4 is recommended for performance reasons.  Configure and install as follows:

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

If the configuration fails consult the **PETSc** user documentation.

Upon completion of the configure step, **PETSc** will provide the make command
that you should use. Execute this command.  After a successful make, **PETSc** should indicate the command used to install the header files and libraries. Execute this command.

To test whether **PETSc** has been installed correctly execute:

    $ make test

Finally, set the PETSC_ROOT environment variables:

    $ export PETSC_ROOT=/Path/to/PETSC

#### Step 5 Install the Visualization Tool Kit

Download and extract **VTK** from https://vtk.org.  **VTK** v8.2.0+ is recommended.  Configure and install as follows:

    $ mkdir vtk-build
    $ cd vtk-build
    $ cmake -DCMAKE_INSTALL_PREFIX=/Path/to/VTK \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DVTK_Group_MPI:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
    /Path/to/VTK/source

    $ make
    $ make install

Set the VTK_DIR environment variables:

    $ export VTK_DIR=/Path/to/VTK

#### Step 6 Configure and Build ChiTech

Clone the **ChiTech** repository:

    $ git clone https://github.com/chi-tech/chi-tech

Run the configure script:

    $ cd chi-tech
    $ ./configure.sh

If all goes well, configure.sh will automatically install the following four packages:
 - readline
 - ncurses
 - Lua 5.3.5
 - Triangle 1.6

If problems are encountered here please see
[troubleshooting installation](CHI_DOC/TroubleShootingInstall.md)

Once configure.sh has completed, build **ChiTech** with:

    $ make

#### Step 7 - Generate the documentation

You can either access the documentation online [here](https://chi-tech.github.io), or generate it locally.

To generate the documentation locally, Google and follow the online instructions for installing
installing doxygen and lua then run:

    $ cd CHI_DOC/
    $ ./YReGenerateDocumentation.sh

Open the documentation index file at CHI_DOC/HTMLdocs/html/index.html.
