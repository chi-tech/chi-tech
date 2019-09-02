#!/bin/sh

# 
# SYNOPSIS:
#   ./configure.sh [clean] [generate] [...]
#
# DESCRIPTION:
#   This script checks the dependencies and installs the required libraries if 
#   they are missing. After that, it generates the CMake build scripts; the 
#   precise behavior of this stage depends on the arguments passed:
#   
#   1. If called without arguments or with a single argument 'generate' (without 
#      quotes), CMake build scripts will be generated in the directory chi_build
#      (which will be created if not present). Additional arguments will be 
#      passed to the CMake generator.
#
#   2. If called with the single argument 'clean', the CMake build directory 
#      ('chi_build') will be removed and the script will stop. If additional 
#      arguments are supplied after 'clean', after removing the build directory, 
#      the script will behave as in case 1 (e.g., to regenerate the CMake build 
#      scripts from scratch, one can use ./configure.sh clean generate).
#

#----- Check if dependencies have been compiled -----
if [ ! -d "CHI_RESOURCES/Dependencies/ncurses" ]; then
  cd "CHI_RESOURCES/Dependencies"
  tar -zxf readline.tar.gz
  tar -zxf ncurses.tar.gz
  tar -zxf lua-5.3.5.tar.gz
  cd readline
  ./configure --prefix=$PWD/build
  make 
  make install
  cd ../ncurses
  ./configure --prefix=$PWD/build
  make 
  make install
  cd ../
  export LIBRARY_PATH=$LIBRARY_PATH:"$PWD/readline/build/lib"
  export LIBRARY_PATH=$LIBRARY_PATH:"$PWD/ncurses/build/lib"
  export CPATH=$CPATH:"$PWD/readline/build/include"

  cd lua-5.3.5
  if [ "$(uname)" == "Darwin" ]; then
      make macosx MYLIBS=-lncurses      
  elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
      make linux MYLIBS=-lncurses
  elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
      make mingw MYLIBS=-lncurses
  elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
      make mingw MYLIBS=-lncurses
  fi
  
  make local

  cd ../
  if [ ! -d "triangle" ]; then
    if [ "$(uname)" == "Darwin" ]; then
        tar -zxf trianglemac.tar.gz  triangle/  
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        tar -zxf triangle.tar.gz
    elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
        tar -zxf triangle.tar.gz
    elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
        tar -zxf triangle.tar.gz
    fi
  fi
  cd triangle
  echo "Building Triangle"
  make triangle
  make trilibrary
  cd ../../../
fi

CMAKE_ARGS=     # additional CMake arguments
DO_CLEAN="No"   # if yes, remove the chi_build directory before generating CMake
DO_CMAKE="Yes"  # by default, run CMake generator after configuring

# Go through all arguments
for arg in "$@"; do

  if [ "$arg" = "clean" ]; then
    # The build directory will be cleaned
    DO_CLEAN="Yes"
    
    if [ $# -eq 1 ]; then
      # If "clean" is the only argument passed
      DO_CMAKE="No"
    fi

  elif [ "$arg" = "generate" ]; then
    DO_CMAKE="Yes"

  else
    # Additional CMake arguments
    CMAKE_ARGS="${CMAKE_ARGS} ${arg}"

  fi

done

if [ $DO_CLEAN = "Yes" ]; then
  # If chi_build directory exists, remove it
  if [ -d "chi_build" ]; then
    rm -r chi_build
  fi
fi

if [ $DO_CMAKE = "Yes" ]; then
  mkdir -p chi_build
  cd chi_build
  cmake "$CMAKE_ARGS" ../
  cd ..
fi
