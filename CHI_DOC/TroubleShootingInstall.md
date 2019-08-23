# Troubleshooting installation

### Compilation cannot find lua headers
The solution to this could be to manually install the lua library.Go to the
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

Compile *lua* (exchange "linux" with "macosx") if you are on a mac

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
    


### If you need to install OpenGL (if graphics needed) (NOT ON MACOSX)

If you want to compile the graphics utilities then make sure the 
*set(CHI_USEGRAPHICS ON)* option is set to *ON* inside the *CMakeLists.txt* file

    sudo apt-get install libx11-dev
    sudo apt-get install libglapi-mesa
    sudo apt-get install mesa-utils
    sudo apt-get install mesa-utils-extra

    sudo apt-get install aptitude

This next command will require a little bit of creativity. The default
Ubuntu 18.04 distribution seems to have been broken. Aptitude will help to prefix
it but the right options have to be chosen

    sudo aptitude install libgl1-mesa-dev

If it gives an option of "leave as current version", select "no". If it gives an option of "downgrade", then
select "yes".

Do the same for the these too (but they shouldn't be broken now):

    sudo aptitude install libglu1-mesa-dev
    sudo aptitude install libglu1-mesa-dev
    sudo aptitude install libglew-dev
