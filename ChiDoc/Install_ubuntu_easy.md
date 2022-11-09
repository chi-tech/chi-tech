# Easy install on Ubuntu Machines

The following instructions were tested on Ubuntu 18.04 LTS and newer LTS (22.04 currently); other Linux
distributions might require some minor tweaking.

Some portions of this document indicate the use of `sudo`. If you do not have 
administrative privileges then have your system administrator assist you in 
these steps.

### Step 1 - Installing GCC, GFortran and the basic environment

- GCC is used to build and install ChiTech.
- GFortran and Python are used during the installation of PETSc
(which ChiTech uses as a linear algebra backend) 
- OpenGL has become optional and VTK can be used, through ViSiT and Paraview, for visualization.

The above packages will therefore need to be installed.

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
updated in which case just run ```sudo apt update```, possibly followed by 
```sudo apt upgrade```.

Now, install the remaining packages needed to build ChiTech and its dependencies:

```bash
sudo apt-get install cmake python3 git zlib1g-dev libx11-dev unzip
```
(note the use of python3)



***TODO: this needs to be fixed @JANV ***
<u>NOTE</u>: *The recommended version of PETSc requires Python 2.6+; Python 3.x, typically
pre-installed on modern Linux systems, will not work (Step 4 below provides more
details if you need to use `python3`) .*

<u>NOTE</u>: If you want to install the *optional* OpenGL package for VTK, do this
```bash
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
```

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

**Important:**  If you want to contribute to **ChiTech**, it is strongly recommended to first fork the **ChiTech** repository into your own Git account and then to clone your fork. 

Clone the **ChiTech** repository.  Go the folder where you want to keep ChiTech relevant stuff:
```bash
    $ git clone https://github.com/chi-tech/chi-tech
```
or
```bash
    $ git clone https://github.com/YOUR-NAME/chi-tech
```

**Important:** We recommend building all the dependencies into a separate folder. For instance,
```bash
    $ mkdir dependencies
```

Before building the dependencies, you need to export a few variables for the PETSc install:
```bash
    $ export CC=mpicc
    $ export CXX=mpicxx
    $ export FC=mpifort
```

Go to the chi-tech folder you have just cloned and type:
```bash
    $ cd chi-tech
    $ python3 ChiResources/configure_dependencies.py ../dependencies
```
The configure script will attempt to download and install all the necessary 
dependencies **and may take a long time**

### Step 4 - Configure environment

The next step in this process is to setup the environment variables for compiling
ChiTech.

```bash
    $. ./chi-dependencies/configure_deproots.sh
```
**Note:** You can replace ```$. ``` in the above with ```$source ```

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
    $ python3 ChiTest/Z_Run_all.py
```

### Step 7 - ChiTech documentation

You can either access the documentation online [here](https://chi-tech.github.io), or generate it locally.

To generate the documentation from your local working copy, first make sure
Doxygen, LaTeX, and lua (currently, lua is at version 5.4) are installed:

```bash
sudo apt-get install doxygen 
sudo apt-get install texlive
sudo apt install texlive-font-utils
sudo apt-get install lua5.4
```

The documentation is contained in the *CHI_DOC* folder and can be generated
using a script provided in that folder:

```bash
cd ChiDoc/
./YReGenerateDocumentation.sh
```

Once finished, you can view the generated documentation by opening

```bash
ChiDoc/HTMLdocs/html/index.html
```

in your favorite browser.