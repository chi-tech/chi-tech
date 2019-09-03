# Troubleshooting installation

### Compilation cannot find lua headers
The solution to this could be to manually install the lua library. Go to the
*Dependencies* folder:

    cd CHI_RESOURCES/Dependencies

Extract *readline*, *ncurses* and *lua-5.3.5*:

    tar -zxvf readline.tar.gz
    tar -zxvf ncurses.tar.gz
    tar -zxvf lua-5.3.5.tar.gz

Compile *readline*

    cd readline
    ./configure --prefix=$PWD/build
    make
    make install

Compile *ncurses*

    cd ../ncurses
    ./configure --prefix=$PWD/build
    make
    make install

Compile *lua* (exchange "linux" with "macosx" if you are on a mac)

    cd ../
    export LIBRARY_PATH=$LIBRARY_PATH:"$PWD/readline/build/lib"
    export LIBRARY_PATH=$LIBRARY_PATH:"$PWD/ncurses/build/lib"
    export CPATH=$CPATH:"$PWD/readline/build/include"
    
    cd lua-5.3.5
    make linux MYLIBS=-lncurses
    make local

### Compilation cannot find Triangle headers or library
Try manually installing it. Go to *Dependencies* folder:

    cd CHI_RESOURCES/Dependencies

Extract the triangle archive

    tar -zxf triangle.tar.gz

Now change directory to this folder and make the library:

    cd triangle/
    make trilibrary