#make clean
#!/bin/sh

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
  if [ "$(uname)" == "Darwin" ]; then
      tar -zxf trianglemac.tar.gz  triangle/  
  elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
      tar -zxf triangle.tar.gz
  elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
      tar -zxf triangle.tar.gz
  elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
      tar -zxf triangle.tar.gz
  fi
  cd triangle
  make
  make trilibrary
  cd ../../../
fi

GR=
DO_CMAKE="Yes"
#===================== If arguments passed
if [ $# -ne 0 ]; then


  if [ "$1" = "G" ]; then
    GR="-DUSE_GRAPHICS=ON"
  fi

  #If two arguments provided
  if [ $# -eq 2 ]; then
    #If 2nd arg is "clean"
    if [ "$2" = "clean" ]; then
      #If directory exists
      if [ -d "chi_build" ]; then
        rm -r chi_build
      fi
    fi
  fi

  #If 1st arg is "clean"
  if [ "$1" = "clean" ]; then
    #If directory exists
    if [ -d "chi_build" ]; then
      rm -r chi_build
    fi
    DO_CMAKE="No"
  fi


fi

if [ "$DO_CMAKE" = "Yes" ]; then
  mkdir -p chi_build
  cd chi_build
  cmake $GR ../
  cd ..
fi


