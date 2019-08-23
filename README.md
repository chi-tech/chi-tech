# Chi-Tech2 #

## Compiling on Personal Linux Machines and MacOS ##

### Step 1 - Installing gcc

Check to see if gcc installed

    gcc --version

This should display something like this:

   $ gcc --version
   gcc (Ubuntu 7.3.0-30ubuntu1~18.04.york0) 7.3.0
   Copyright (C) 2017 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If this is not displayed google how to do it and get it done before proceeding.

### Step 2 - An MPI flavor

Download and install either OpenMPI or MPICH. If you have MOOSE or deal.ii 
installed then you are probably already set to go.

### Step 3 - Boost 1.63+
Download and install boost. Eventually ChiTech needs the location where the 
"include" and "lib" directories are so just follow online instructions for this.

### Step 4 - PetSc

The best performance thus far tested is with 
[petsc version 3.9.4](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz)
it is recommended to use this version.

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

### Step 5 - Install cmake and X11

    sudo apt-get install cmake libx11-dev
    
On MacOs you will have to find the equivalent using homebrew. 
This can easily be googled.

### Step 6 Configure ChiTech
Go to the chi_tech_2 folder and type

    ./configure.sh

If all goes well it will automatically install 4 things
 - The readline library,
 - The ncurses library,
 - Lua 5.3.5, and
 - Triangle 1.6
 
If problems are encountered here please see 
[troubleshooting installation](CHI_DOC/TroubleShootingInstall.md)



### Step 7 - Build ChiTech

Navigate to the main directory (i.e. *chi_tech_2/*)

    ./configure.sh
    make -j4

You can also use -j8 even if you don't have 8 processors, the make command 
will use threading where possible.

### Step 8 - Generate the documentation

Google and follow online instructions for installing doxygen then

    cd CHI_DOC/
    ./YReGenerateDocumentation.sh 
    
Open the documentation index file at 

    CHI_DOC/HTMLdocs/html/index.html