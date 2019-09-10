## Compiling on macOS
___
### Compiling with the Apple Tdoolchain
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

    $ make all test

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

Set the VTK_INSTALL environment variables:

    $ export VTK_INSTALL=/Path/to/VTK

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

Google and follow online instructions for installing doxygen and lua then run:

    $ cd CHI_DOC/
    $ ./YReGenerateDocumentation.sh

Open the documentation index file at CHI_DOC/HTMLdocs/html/index.html.
___
### Compiling with the GNU toolchain
___
#### Step 1 - Installing gcc, gfortran and the basic environment

GCC is used to install and run ChiTech. Gfortran is used
during the installation of PETSc and therefore also needs to
be installed.

Check to see if gcc is installed

    gcc --version

This should display something like this:

    $ gcc --version
    gcc (Ubuntu 7.3.0-30ubuntu1~18.04.york0) 7.3.0
    Copyright (C) 2017 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If this is not displayed then install gcc using the following command:

    sudo apt-get install build-essential gfortran

If you still do not get the appropriate version when running either ``gcc --version`` or ``gfortran --version`` then there are many online resources which may be able to assist you.

Install other packages:

    sudo apt-get install cmake \
    libx11-dev \
    doxygen \
    zlib1g-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    python python3 \
    git


#### Step 2 - An MPI flavor

Install either OpenMPI or MPICH. If you have MOOSE or deal.ii
installed then you are probably already set to go. **MPICH** is recommended for
better performance.

    sudo apt-get install mpich

To check if this is working properly, simply do the mpi version of checking for
a c++ compiler:

    mpicc --version

Which should display the same message the gcc call did, i.e.

    $ mpicc --version
    gcc (Ubuntu 7.3.0-30ubuntu1~18.04.york0) 7.3.0
    Copyright (C) 2017 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#### Step 3 - Boost 1.63+
Download and unpack boost. Eventually ChiTech needs the location where the
"include" directory is so just follow online instructions for this.

Alternatively, just make a *projects* directory and download boost there,
like the following:

    mkdir projects
    cd projects
    mkdir boost
    cd boost
    wget https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.gz
    tar -zxf boost_1_71_0.tar.gz
    cd boost_1_17_0
    ./bootstrap.sh
    ./b2 headers
    mkdir include
    cp -r boost include/

The current directory will now be your *BOOST_ROOT* folder.

    export BOOST_ROOT=$PWD

You probably want to permanently add this to your bash profile file.
This is done in different ways depending on your platform.

#### Step 4 - PetSc

The best performance thus far tested is with
[petsc version 3.9.4](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz)
it is recommended to use this version.

Return to your *projects* folder (or whatever you chose to place stuff). Run
the following

    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz
    tar -zxf petsc-3.9.4.tar.gz
    cd petsc-3.9.4

Download and extract the archive file to a folder of your choice then navigate
into the directory containing the "configure" script and execute the following:

    ./configure --prefix=$PWD/install \
    --download-hypre=1 \
    --with-ssl=0 \
    --with-debugging=0 \
    --with-pic=1 \
    --with-shared-libraries=1 \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-parmetis=1 \
    --download-superlu_dist=1 \
    --with-cxx-dialect=C++11 \
    CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 \
    CFLAGS='-fPIC -fopenmp' \
    CXXFLAGS='-fPIC -fopenmp' \
    FFLAGS='-fPIC -fopenmp' \
    FCFLAGS='-fPIC -fopenmp' \
    F90FLAGS='-fPIC -fopenmp' \
    F77FLAGS='-fPIC -fopenmp' \
    COPTFLAGS='-O3 -march=native -mtune=native' \
    CXXOPTFLAGS='-O3 -march=native -mtune=native' \
    FOPTFLAGS='-O3 -march=native -mtune=native' \
    PETSC_DIR=$PWD

If the configuration fails then consult PetSc's user documentation.

Upon completion of the configure step, petsc will provide the make command
that you should use. Execute this command.

After a successful make petsc should indicate the command used to install
the header files and libraries. Execute the associated command as well.

To test whether the system has been installed correctly execute:

    make all test

The final step is to generate the environment variable for the install
directory of PETSc:

    export PETSC_ROOT=$PWD/install

Again, this is also something you'd like to add to your bash profile.

#### Step 5 Install the Visualization Tool Kit

In your projects folder install VTK using the following commands:

    mkdir VTK
    cd VTK
    wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
    tar -zxf VTK-8.2.0.tar.gz
    cd VTK-8.2.0
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DVTK_Group_MPI:BOOL=ON \
    -DCMAKE_BUILD_TYPE=Release \
     ../

This will take a while after which you need to execute:

    make -j4 && make install

This will also take a while. Finally export the *VTK_ROOT* environment
variable as well as the version:

    export VTK_ROOT=$PWD/../install
    export VTK_VERSION=-8.2

Since we will be dependent on shared libraries from VTK we will also export this
to our shared library path

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VTK_ROOT/lib

#### Step 6 Configure ChiTech
Clone the repository. Go the folder where you want to keep ChiTech relevant stuff:

    git clone https://github.com/chi-tech/chi-tech
    cd chi-tech

Go to the chi_tech_2 folder and type

    ./configure.sh

If all goes well it will automatically install 4 things
 - The readline library,
 - The ncurses library,
 - Lua 5.3.5, and
 - Triangle 1.6

If problems are encountered here please see
[troubleshooting installation](CHI_DOC/TroubleShootingInstall.md)

#### Step 7 - Build ChiTech

Navigate to the main directory (i.e. *chi_tech_2/*)

    ./configure.sh
    make -j4

You can also use -j8 even if you don't have 8 processors, the make command
will use threading where possible.

#### Step 8 - Generate the documentation

Google and follow online instructions for installing doxygen and lua then run

    cd CHI_DOC/
    ./YReGenerateDocumentation.sh

Open the documentation index file at

    CHI_DOC/HTMLdocs/html/index.html