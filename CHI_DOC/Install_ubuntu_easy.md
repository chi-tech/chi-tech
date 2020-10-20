# Easy install on Ubuntu Machines

The following instructions were tested on Ubuntu 18.04 LTS; other Linux
distributions might require some minor tweaking.

Some portions of this document indicate the use of `sudo`. If you do not have 
administrative privileges then have you system administrator assist you in 
these steps.

### Step 1 - Installing GCC, GFortran and the basic environment

GCC is used to build and install ChiTech.
GFortran and Python are used during the installation of PETSc
(which ChiTech uses as a linear algebra backend) and
OpenGL is required by VTK (used by ChiTech for visualization).
These packages will therefore also need to be installed.

Check to see if gcc is installed

```bash
gcc --version
```

This should display something like this:

    $ gcc --version
    gcc (Ubuntu 7.3.0-30ubuntu1~18.04.york0) 7.3.0
    Copyright (C) 2017 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If this is not displayed then install gcc using the following command:

```bash
sudo apt-get install build-essential gfortran
```

If you still do not get the appropriate version when running either ``gcc --version``
or ``gfortran --version`` then there are many online resources which may be able to
assist you. One issue commonly seen is that the repositories on your system are not
updated in which case just run ```sudo apt update``` followed by 
```sudo apt upgrade```.

Now, install the remaining packages needed to build ChiTech and its dependencies:

```bash
sudo apt-get install cmake python git zlib1g-dev libx11-dev unzip
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
```

<u>NOTE</u>: *The recommended version of PETSc requires Python 2.6+; Python 3.x, typically
pre-installed on modern Linux systems, will not work (Step 4 below provides more
details if you need to use `python3`) .*

<u>NOTE</u>: *The second line is to install OpenGL for VTK.*

### Step 2 - An MPI flavor

Install either OpenMPI or MPICH. If you have MOOSE or deal.ii installed then you
are probably already set to go. **MPICH** is recommended for better performance.

```bash
sudo apt-get install mpich
```

To check if this is working properly, simply do the mpi version of checking for
a C++ compiler:

```bash
mpicc --version
```

Which should display the same message the gcc call did, i.e.,

    $ mpicc --version
    gcc (Ubuntu 7.3.0-30ubuntu1~18.04.york0) 7.3.0
    Copyright (C) 2017 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

### Step 3 - Clone ChiTech

Clone the **ChiTech** repository.  Go the folder where you want to keep ChiTech relevant stuff:
```bash
    $ git clone https://github.com/chi-tech/chi-tech
```

Go to the chi-tech folder and type:
```bash
    $ cd chi-tech
    $ python3 CHI_RESOURCES/configure_dependencies.py
```
The configure script will attempt to download and install all the necessary 
dependencies **and may take a long time**

### Step 4 - Configure environment

The next step in this process is to setup the environment variables for compiling
ChiTech.

```bash
    $. ./chi-dependencies/configure_deproots.sh
```


### Step 5 - Build ChiTech

To compile ChiTech now just execute:
```bash
    $ ./configure.sh
```
The configure script will generate the CMake build scripts.

In the main directory (i.e. *chi-tech/*), execute:
```bash
    $ make -j4
```

You can also use -j8 even if you don't have 8 processors, the make command
will use threading where possible.

Whenever recompilation is needed or the configuration has changed,
you can run
```bash
    $ ./configure.sh clean
    $ ./configure.sh
    $ make -j4
```

### Step 6 - Run regression tests

To check if the code compiled correctly execute the test scripts:

```bash
    $ python3 CHI_TEST/Z_Run_all.py
```

### Step 7 - ChiTech documentation

You can either access the documentation online [here](https://chi-tech.github.io), or generate it locally.

To generate the documentation from your local working copy, first make sure
Doxygen and LaTeX are installed:

```bash
sudo apt-get install doxygen texlive
```

The documentation is contained in the *CHI_DOC* folder and can be generated
using a script provided in that folder:

```bash
cd CHI_DOC/
./YReGenerateDocumentation.sh
```

Once finished, you can view the generated documentation by opening

```bash
CHI_DOC/HTMLdocs/html/index.html
```

in your favorite browser.